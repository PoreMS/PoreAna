################################################################################
# Sample                                                                       #
#                                                                              #
"""Sample functions to be run on clusters."""
################################################################################

import math
import multiprocessing as mp
import sys

import chemfiles as cf
import numpy as np
import porems as pms
import scipy

import poreana.geometry as geometry
import poreana.utils as utils


class Sample:
    """This class samples a trajectory to determine different properties.
    Different properties can be initialized to be run at the same time during
    the sampling run. The output is stored in form of pickle files for later
    calculation using methods provided in the package.

    It is advisable to run the sampling on a cluster due to a high time and
    resource consumption.

    The system can either be a pore system - variable **system** is a file link
    to the *pore_system* object file - or a simple simulation box - variable
    **system** is a list containing the dimensions in nano meter.

    Parameters
    ----------
    system : string, list
        Link to poresystem object file or a list of box dimensions (nm)
    link_traj : string
        Link to trajectory file (trr or xtc)
    mol : Molecule
        Molecule to calculate the density for
    atoms : list, optional
        List of atom names, leave empty for whole molecule
    masses : list, optional
        List of atom masses, leave empty to read molecule object masses
    entry : float, optional
        Remove pore entrance from calculation
    frame_end : int, optional
        Set an end frame to stop the analysis earlier
    """

    def __init__(
        self, system, link_traj, mol, atoms=[], masses=[], entry=0.5, frame_end=-1
    ):
        self._pore = (
            utils.load(system, file_type="yml") if isinstance(system, str) else None
        )
        self._box = system if isinstance(system, list) else []
        self._traj = link_traj
        self._mol = mol
        self._atoms = atoms
        self._masses = masses
        self._entry = entry

        self._is_density = False
        self._is_gyration = False
        self._is_angle = False
        self._is_diffusion_bin = False
        self._is_diffusion_mc = False
        self._is_diffusion_vacf = False
        self._is_numpy = False

        # Resolve atom indices from names
        self._atoms = (
            [atom.get_name() for atom in mol.get_atom_list()]
            if not self._atoms
            else self._atoms
        )
        self._atoms = [
            i
            for i in range(mol.get_num())
            if mol.get_atom_list()[i].get_name() in self._atoms
        ]

        # Resolve masses
        if not self._masses:
            if len(self._atoms) == mol.get_num():
                self._masses = mol.get_masses()
            elif len(self._atoms) == 1:
                self._masses = [1]
        self._sum_masses = sum(self._masses)

        if self._atoms and len(self._masses) != len(self._atoms):
            print("Length of variables *atoms* and *masses* do not match!")
            return

        # Pre-compute numpy arrays for hot-loop efficiency
        self._masses_arr = np.array(self._masses, dtype=float)
        self._atoms_per_mol = mol.get_num()

        # Get trajectory metadata
        traj = cf.Trajectory(self._traj)
        self._num_frame = traj.nsteps if frame_end == -1 else frame_end

        frame = traj.read()
        num_res = len(frame.topology.atoms) / mol.get_num()

        if abs(int(num_res) - num_res) >= 1e-5:
            print("Number of atoms is inconsistent with number of residues.")
            return

        # Build residue → atom-index mapping as numpy arrays for fast slicing
        n_res = int(num_res)
        self._res_list = {
            res_id: np.array(
                [
                    res_id * mol.get_num() + a
                    for a in range(mol.get_num())
                    if a in self._atoms
                ]
            )
            for res_id in range(n_res)
        }

        # Extract pore shape properties
        self._pore_props = {}
        if self._pore:
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    p = self._pore[pore_id]
                    self._pore_props[pore_id] = {
                        "type": p["shape"],
                        "focal": p["parameter"]["centroid"],
                        "length": p["parameter"]["length"],
                        "diam": p["diameter"],
                    }
            self._pore_props["box"] = {
                "dimensions": self._pore["system"]["dimensions"],
                "res": self._pore["system"]["reservoir"],
            }

    ########
    # Bins #
    ########
    def _bin_in(self, bin_num):
        """Radial bin structure for the pore interior."""
        width = {}
        for pore_id in self._pore_props:
            if pore_id[:5] == "shape":
                r_max = self._pore_props[pore_id]["diam"] / 2
                width[pore_id] = [r_max / bin_num * x for x in range(bin_num + 2)]
        return {"width": width, "bins": [0] * (bin_num + 1)}

    def _bin_in_slit(self, bin_num):
        """Bin structure across the slit pore height."""
        y_max = self._pore_props["box"]["dimensions"][1]
        width = {}
        for pore_id in self._pore_props:
            if pore_id[:5] == "shape":
                width[pore_id] = [y_max / bin_num * x for x in range(bin_num + 2)]
        return {"width": width, "bins": [0] * (bin_num + 1)}

    def _bin_in_const_A(self, bin_num):
        """Radial bin structure with equal bin area for the pore interior."""
        width = {}
        bins = None
        for pore_id in self._pore:
            if pore_id[:5] == "shape":
                diam = self._pore_props[pore_id]["diam"]
                n = bin_num + 2
                pore_surf = diam**2
                surf_per_bin = pore_surf / n

                matrix_bins = np.zeros((n, n))
                for i in range(n):
                    if i == 0:
                        matrix_bins[i, 0] = 1
                    elif i == n - 1:
                        matrix_bins[i, -2] = 1
                    else:
                        matrix_bins[i, i] = 1
                        matrix_bins[i, i - 1] = -1

                res_vec = [surf_per_bin] * n
                res_vec[-1] = -surf_per_bin + (diam / 2) ** 2

                x = scipy.sparse.linalg.lsmr(matrix_bins, res_vec)[0]
                x[-1] = (diam / 2) ** 2
                r_vals = list(np.sqrt(x)[:-1])
                r_vals.insert(0, 0)
                width[pore_id] = r_vals
                bins = [0] * (n + 1)

        return {"width": width, "bins": bins}

    def _bin_ex(self, bin_num):
        """Bin structure for the pore exterior / reservoir."""
        length = (
            self._pore_props["box"]["res"]
            if self._pore_props
            else self._box[self._dens_inp["direction"]]
        )
        width = list(np.linspace(0, length, bin_num + 1))
        return {"width": width, "bins": [0] * (bin_num + 1)}

    def _bin_window(self, bin_num, len_window):
        """Window bin structure for bin-diffusion inside the pore."""
        width = {}
        for pore_id in self._pore_props:
            if pore_id[:5] == "shape":
                r_max = self._pore_props[pore_id]["diam"] / 2
                width[pore_id] = [r_max / bin_num * x for x in range(bin_num + 2)]
        return {"width": width, "bins": [[0] * len_window for _ in range(bin_num + 1)]}

    def _bin_mc(self, bin_num, direction):
        """Bin structure spanning the full simulation length for MC diffusion."""
        z_length = (
            self._pore_props["box"]["dimensions"][direction]
            if self._pore
            else self._box[direction]
        )
        bins = list(np.linspace(0, z_length, bin_num + 1))
        return {"bins": bins}

    def _bin_pore(self, bin_num):
        """Radial bin structure for the pore interior based on pore diameter."""
        bins = {}
        for pore_id in self._pore:
            if pore_id[:5] == "shape":
                ptype = self._pore_props[pore_id]["type"]
                if ptype == "CYLINDER":
                    bins[pore_id] = np.linspace(
                        0, self._pore_props[pore_id]["diam"] / 2, bin_num + 1
                    )
                else:
                    return f"{ptype} radial binning not implemented yet."
        return {"bins": bins}

    ###########
    # Density #
    ###########
    def init_density(
        self,
        link_out,
        bin_num=150,
        remove_pore_from_res=False,
        bin_const_A=False,
        avg_slit=True,
        direction=2,
    ):
        """Enable density sampling routine.

        Parameters
        ----------
        link_out : string
            Link to output hdf5, obj or yml data file
        bin_num : integer, optional
            Number of bins to be used
        remove_pore_from_res : bool, optional
            True to remove an extended pore volume from the reservoirs
        bin_const_A : bool, optional
            If true, all radial bins will have the same surface area
        avg_slit : bool, optional
            If False the density profile over the slit pore height is calculated
        direction : int, optional
            Direction for box-system density (x=0, y=1, z=2)
        """
        self._is_density = True
        self._dens_inp = {
            "output": link_out,
            "bin_num": bin_num,
            "remove_pore_from_res": remove_pore_from_res,
            "bin_const_A": bin_const_A,
            "avg_slit": avg_slit,
            "direction": direction,
        }

    def _density_data(self):
        """Create density data structure."""
        bin_num = self._dens_inp["bin_num"]
        ex = self._bin_ex(bin_num)
        data = {"ex_width": ex["width"], "ex": ex["bins"]}

        if self._pore:
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    data[pore_id] = {}
                    if self._dens_inp["bin_const_A"]:
                        inner = self._bin_in_const_A(bin_num)
                        data[pore_id]["in_width"] = inner["width"][pore_id]
                        data[pore_id]["in"] = inner["bins"]
                    elif not self._dens_inp["avg_slit"]:
                        inner = self._bin_in_slit(bin_num)
                        data[pore_id]["in_width"] = inner["width"][pore_id]
                        data[pore_id]["in"] = inner["bins"]
                    else:
                        inner = self._bin_in(bin_num)
                        data[pore_id]["in_width"] = inner["width"][pore_id]
                        data[pore_id]["in"] = inner["bins"]
        return data

    def _density(self, data, region, dist, com, pore_id):
        """Sample density inside and outside of the pore."""
        bin_num = self._dens_inp["bin_num"]

        if region == "in" and pore_id != 0:
            if self._dens_inp["bin_const_A"]:
                index = np.digitize(dist[pore_id], data[pore_id]["in_width"][1:])
            elif not self._dens_inp["avg_slit"]:
                index = np.digitize(com[1], data[pore_id]["in_width"][1:])
            else:
                index = int(dist[pore_id] / data[pore_id]["in_width"][1])

            if index <= bin_num:
                data[pore_id]["in"][index] += 1

        elif region == "ex":
            length = (
                abs(com[2] - self._pore_props["box"]["dimensions"][2])
                if self._pore and com[2] > self._pore_props["box"]["dimensions"][2] / 2
                else com[self._dens_inp["direction"]]
            )
            index = int(length / data["ex_width"][1])

            if (
                self._pore
                and self._dens_inp["remove_pore_from_res"]
                and pore_id == "shape_00"
            ):
                is_add = (
                    index <= bin_num
                    and dist[pore_id] > self._pore_props[pore_id]["diam"] / 2
                    and com[2] <= self._pore_props["box"]["dimensions"][2]
                )
            else:
                is_add = index <= bin_num

            if is_add:
                data["ex"][index] += 1

    ############
    # Gyration #
    ############
    def init_gyration(self, link_out, bin_num=150):
        """Enable gyration sampling routine.

        Parameters
        ----------
        link_out : string
            Link to output hdf5, obj or yml data file
        bin_num : integer, optional
            Number of bins to be used
        """
        self._is_gyration = True
        self._gyr_inp = {"output": link_out, "bin_num": bin_num}

    def _gyration_data(self):
        """Create gyration data structure."""
        bin_num = self._gyr_inp["bin_num"]
        ex = self._bin_ex(bin_num)
        data = {"ex_width": ex["width"], "ex": ex["bins"]}

        if self._pore:
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    inner = self._bin_in(bin_num)
                    data[pore_id] = {
                        "in_width": inner["width"][pore_id],
                        "in": inner["bins"],
                    }
        return data

    def _gyration(self, data, region, dist, com, pos, pore_id):
        r"""Sample the radius of gyration inside and outside of the pore.

        .. math::

            R_g=\left(\frac{\sum_i\|\boldsymbol{r}_i\|^2m_i}{\sum_im_i}\right)^{\frac{1}{2}}

        Parameters
        ----------
        data : dictionary
            Data dictionary with bins for pore interior and exterior
        region : string
            'in' or 'ex'
        dist : dict
            Distance of CoM to pore axis per pore shape
        com : array-like
            Centre of mass of the current molecule
        pos : ndarray
            Atom positions of the current molecule, shape (n_atoms, 3)
        pore_id : str or int
            Active pore id, or 0 for exterior
        """
        bin_num = self._gyr_inp["bin_num"]

        diffs = pos - np.asarray(com)
        r_sq = np.einsum("ij,ij->i", diffs, diffs)
        r_g = float(np.sqrt(np.dot(r_sq, self._masses_arr) / self._sum_masses))

        if region == "in" and pore_id != 0:
            index = int(dist[pore_id] / data[pore_id]["in_width"][1])
            if index <= bin_num:
                data[pore_id]["in"][index] += r_g

        elif region == "ex":
            length = (
                abs(com[2] - self._pore_props["box"]["dimensions"][2])
                if self._pore and com[2] > self._pore_props["box"]["dimensions"][2] / 2
                else com[2]
            )
            index = int(length / data["ex_width"][1])

            if self._pore and pore_id != 0:
                is_add = (
                    index <= bin_num
                    and dist[pore_id] > self._pore_props[pore_id]["diam"] / 2
                    and com[2] <= self._pore_props["box"]["dimensions"][2]
                )
            else:
                is_add = index <= bin_num

            if is_add:
                data["ex"][index] += r_g

    #########
    # Angle #
    #########
    def init_angle(self, link_out, vector_atoms, bin_num=150, normals={}):
        """Enable angle sampling routine.

        Parameters
        ----------
        link_out : string
            Link to output hdf5, obj or yml data file
        vector_atoms : list
            List of two atom ids defining the molecule vector
        bin_num : integer, optional
            Number of bins to be used
        normals : dictionary, optional
            Surface normal vector functions for interior *in* and exterior *ex*
        """
        self._is_angle = True

        if not normals:
            normals = {}
            if self._pore:
                for pore_id in self._pore:
                    if pore_id[:5] == "shape":
                        normals[pore_id] = {}
                        if self._pore_props[pore_id]["type"] == "CYLINDER":
                            shape = pms.Cylinder(
                                {
                                    "centroid": self._pore_props[pore_id]["focal"],
                                    "central": [0, 0, 1],
                                    "length": self._pore_props["box"]["dimensions"][2],
                                    "diameter": self._pore_props[pore_id]["diam"],
                                }
                            )

                            def normal_in(pos, _s=shape):
                                return _s.normal(pos)

                            def normal_ex(pos, _pp=self._pore_props):
                                return (
                                    [0, 0, -1]
                                    if pos[2]
                                    < (_pp["box"]["dimensions"][2] - _pp["box"]["res"])
                                    else [0, 0, 1]
                                )

                            normals[pore_id] = {"in": normal_in}
                            normals["ex"] = normal_ex
                        else:
                            print(
                                "Angle: Shape normal not predefined yet. Please set the 'normals' variable..."
                            )
                            return
            else:
                normals = {
                    "in": lambda pos: [0, 0, 1],
                    "ex": lambda pos: [0, 0, 1],
                }
        self._angle_normals = normals
        self._angle_inp = {
            "output": link_out,
            "vector_atoms": vector_atoms,
            "bin_num": bin_num,
        }

    def _angle_data(self):
        """Create angle data structure."""
        bin_num = self._angle_inp["bin_num"]
        ex = self._bin_ex(bin_num)
        data = {"ex_width": ex["width"], "ex": ex["bins"]}

        if self._pore:
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    inner = self._bin_in(bin_num)
                    data[pore_id] = {
                        "in_width": inner["width"][pore_id],
                        "in": inner["bins"],
                    }
        return data

    def _angle(self, data, region, dist, com, pos, pore_id):
        """Sample the molecule-surface angle inside and outside the pore.

        Parameters
        ----------
        data : dictionary
            Data dictionary with bins
        region : string
            'in' or 'ex'
        dist : dict
            Distance of CoM to pore axis per shape
        com : array-like
            Centre of mass of the current molecule
        pos : ndarray
            Atom positions of the current molecule, shape (n_atoms, 3)
        pore_id : str or int
            Active pore id, or 0 for exterior
        """
        bin_num = self._angle_inp["bin_num"]
        vector_atoms = self._angle_inp["vector_atoms"]
        normals = self._angle_normals

        if region == "in" and pore_id != 0:
            vec = pos[vector_atoms[1]] - pos[vector_atoms[0]]
            angle_val = geometry.angle(vec, normals[pore_id]["in"](com))
            index = int(dist[pore_id] / data[pore_id]["in_width"][1])
            if index <= bin_num:
                data[pore_id]["in"][index] += angle_val

        elif region == "ex":
            vec = pos[vector_atoms[1]] - pos[vector_atoms[0]]
            angle_val = geometry.angle(vec, normals["ex"](com))
            length = (
                abs(com[2] - self._pore_props["box"]["dimensions"][2])
                if self._pore and com[2] > self._pore_props["box"]["dimensions"][2] / 2
                else com[2]
            )
            index = int(length / data["ex_width"][1])

            if self._pore and pore_id != 0:
                is_add = (
                    index <= bin_num
                    and dist[pore_id] > self._pore_props[pore_id]["diam"] / 2
                    and com[2] <= self._pore_props["box"]["dimensions"][2]
                )
            else:
                is_add = index <= bin_num

            if is_add:
                data["ex"][index] += angle_val

    #############
    # Diffusion #
    #############
    def init_diffusion_bin(
        self,
        link_out,
        bin_num=50,
        len_obs=16e-12,
        len_frame=2e-12,
        len_step=2,
        bin_step_size=1,
    ):
        """Enable bin-diffusion sampling routine.

        Parameters
        ----------
        link_out : string
            Link to output obj or yml data file
        bin_num : integer, optional
            Number of radial bins
        len_obs : float, optional
            Observation window length in seconds
        len_frame : float, optional
            Frame length in seconds
        len_step : integer, optional
            Step size between frames
        bin_step_size : integer, optional
            Allowed radial bin displacement per window
        """
        if self._is_diffusion_mc:
            print("Binning and MC-approaches cannot be run in parallel.")
            return
        if self._is_diffusion_vacf:
            print("Binning and VACF-approaches cannot be run in parallel.")
            return
        if self._is_numpy:
            print("Binning and numpy-approaches cannot be run in parallel.")
            return
        if not self._pore:
            print("Bin diffusion currently only usable for pore system.")
            return
        self._is_diffusion_bin = True

        len_window = len_obs / len_step / len_frame + 1
        if len_window != int(len_window):
            obs_u = (math.ceil(len_window) - 1) * len_step * len_frame
            obs_d = (math.floor(len_window) - 1) * len_step * len_frame
            print(
                f"Observation length not possible. Use len_obs={'%.1e' % obs_u} or {'%.1e' % obs_d}."
            )
            return
        len_window = int(len_window)

        self._diff_bin_inp = {
            "output": link_out,
            "bin_step_size": bin_step_size,
            "bin_num": bin_num,
            "len_step": len_step,
            "len_frame": len_frame,
            "len_window": len_window,
        }

    def _diffusion_bin_data(self):
        """Create bin-diffusion data structure."""
        bin_num = self._diff_bin_inp["bin_num"]
        len_window = self._diff_bin_inp["len_window"]
        data = {}
        for pore_id in self._pore:
            if pore_id[:5] == "shape":
                w = self._bin_window(bin_num, len_window)
                data[pore_id] = {"width": w["width"][pore_id]}
                for key in ("z", "r", "n", "z_tot", "r_tot", "n_tot"):
                    data[pore_id][key] = [[0] * len_window for _ in range(bin_num + 1)]
        return data

    def _diffusion_bin_step(self, idx):
        """List of bin indices within the allowed displacement of idx."""
        s = self._diff_bin_inp["bin_step_size"]
        return list(range(idx + s, idx - 1, -1)) + list(range(idx - 1, idx - s - 1, -1))

    def _diffusion_bin(
        self, data, region, pore_in, dist, com_list, idx_list, res_id, com
    ):
        r"""Sample mean-square displacement for bin-diffusion inside the pore.

        Parameters
        ----------
        data : dictionary
            Data dictionary with axial/radial msd bins
        region : string
            'in' or 'ex'
        pore_in : str or int
            Active pore id
        dist : dict
            Distance of CoM to pore axis per shape
        com_list : dict
            Sliding window of CoM per residue per pore
        idx_list : dict
            Sliding window of bin indices per residue per pore
        res_id : integer
            Current residue id
        com : array-like
            Centre of mass of the current molecule
        """
        bin_num = self._diff_bin_inp["bin_num"]
        len_step = self._diff_bin_inp["len_step"]
        len_window = self._diff_bin_inp["len_window"]

        if region != "in":
            return

        index = math.floor(dist[pore_in] / data[pore_in]["width"][1])
        com_list[pore_in][-1][res_id] = com
        idx_list[pore_in][-1][res_id] = index

        if (
            len(com_list[pore_in]) == len_window * len_step
            and res_id in com_list[pore_in][0]
        ):
            pos_ref = com_list[pore_in][0][res_id]
            idx_ref = idx_list[pore_in][0][res_id]

            msd_z = [0] * len_window
            msd_r = [0] * len_window
            norm = [0] * len_window
            len_msd = 0

            for step in range(0, len_window * len_step, len_step):
                if res_id not in com_list[pore_in][step]:
                    break
                pos_step = com_list[pore_in][step][res_id]
                idx_step = idx_list[pore_in][step][res_id]
                win_idx = step // len_step

                msd_z[win_idx] += (pos_ref[2] - pos_step[2]) ** 2
                msd_r[win_idx] += (
                    geometry.length(
                        geometry.vector(pos_ref, [pos_step[0], pos_step[1], pos_ref[2]])
                    )
                    ** 2
                )
                norm[win_idx] += 1

                if idx_step in self._diffusion_bin_step(idx_ref):
                    len_msd += 1
                else:
                    break

            if idx_ref <= bin_num:
                for i in range(len_window):
                    data[pore_in]["z_tot"][idx_ref][i] += msd_z[i]
                    data[pore_in]["r_tot"][idx_ref][i] += msd_r[i]
                    data[pore_in]["n_tot"][idx_ref][i] += norm[i]
                    if len_msd == len_window:
                        data[pore_in]["z"][idx_ref][i] += msd_z[i]
                        data[pore_in]["r"][idx_ref][i] += msd_r[i]
                        data[pore_in]["n"][idx_ref][i] += norm[i]

    ################
    # MC Diffusion #
    ################
    def init_diffusion_mc(
        self, link_out, len_step, bin_num=100, len_frame=2e-12, direction=2
    ):
        r"""Enable MC-diffusion transition-matrix sampling.

        Parameters
        ----------
        link_out : string
            Link to output obj or yml data file
        len_step : list
            Step lengths (lag times) to sample
        bin_num : integer, optional
            Number of spatial bins
        len_frame : float, optional
            Frame length in seconds
        direction : integer, optional
            Spatial direction (0=x, 1=y, 2=z)
        """
        if self._is_diffusion_bin:
            print("Binning and MC-approaches cannot be run in parallel.")
            return
        if self._is_diffusion_vacf:
            print("VACF and MC-approaches cannot be run in parallel.")
            return
        if self._is_numpy:
            print("Numpy and MC-approaches cannot be run in parallel.")
            return
        if direction not in [0, 1, 2]:
            print(
                "Wrong directional input. Possible inputs are 0 (x-axis), 1 (y-axis), and 2 (z-axis)..."
            )
            return

        self._is_diffusion_mc = True
        bins = self._bin_mc(bin_num, direction)["bins"]
        len_step = sorted(len_step)

        self._diff_mc_inp = {
            "output": link_out,
            "bins": bins,
            "bin_num": bin_num,
            "len_step": len_step,
            "len_frame": len_frame,
            "is_pbc": True,
            "direction": int(direction),
        }

    def _diffusion_mc_data(self):
        """Create MC-diffusion data structure."""
        bin_num = self._diff_mc_inp["bin_num"]
        return {
            step: np.zeros((bin_num + 2, bin_num + 2), int)
            for step in self._diff_mc_inp["len_step"]
        }

    def _diffusion_mc(self, data, idx_list, com, res_id, frame_list, frame_id):
        """Sample the MC-diffusion transition matrix.

        Parameters
        ----------
        data : dictionary
            Transition matrices keyed by step length
        idx_list : list
            Sliding list of bin-index dicts per frame
        com : array-like
            CoM of the current molecule
        res_id : integer
            Current residue id
        frame_list : list
            Frame ids assigned to this worker
        frame_id : integer
            Current frame id
        """
        len_step = self._diff_mc_inp["len_step"]
        bins = self._diff_mc_inp["bins"]
        direction = self._diff_mc_inp["direction"]

        idx_list[-1][res_id] = np.digitize(com[direction], bins)

        if frame_list[0] == 0:
            for step in len_step:
                if len(idx_list) >= step + 1:
                    data[step][idx_list[-1][res_id], idx_list[-(step + 1)][res_id]] += 1

        elif frame_id >= frame_list[0] + self._diff_mc_inp["len_step"][-1]:
            for step in len_step:
                if len(idx_list) >= step + 1:
                    data[step][idx_list[-1][res_id], idx_list[-(step + 1)][res_id]] += 1

    ######################
    # VACF Diffusion     #
    ######################
    def init_diffusion_vacf(
        self,
        link_out: str,
        len_correration=2e-11,
        new_time_origin=2e-13,
        sample_step=20,
        len_frame=1e-15,
        bin_num=32,
        direction=2,
        sample_each_residue=False,
    ):
        """Enable VACF diffusion sampling for local diffusion coefficients.

        Samples the local velocity autocorrelation function (VACF) to compute
        direction-resolved diffusion coefficients per spatial bin. For cylindrical
        pores, ``direction='radial'`` bins in the radial coordinate and computes
        diffusion in cylindrical coordinates (radial, tangential, axial).

        Parameters
        ----------
        link_out : string
            Output file path
        len_correration : float, optional
            Correlation time in seconds
        new_time_origin : float, optional
            Time between successive time origins in seconds
        sample_step : integer, optional
            Frames between velocity samples
        len_frame : float, optional
            Frame length in seconds
        bin_num : integer, optional
            Number of spatial bins
        direction : integer or ``'radial'``, optional
            Binning direction: 0 (x), 1 (y), 2 (z), or ``'radial'`` for
            cylindrical pore systems
        sample_each_residue : bool, optional
            Sample VACF per residue instead of per bin average
        """
        corr_steps = int(round(len_correration / len_frame / sample_step))
        new_time_origin_steps = int(round(new_time_origin / len_frame / sample_step))
        if corr_steps < 1 or new_time_origin_steps < 1:
            print(
                "VACF needs a correlation time longer than one frame. "
                "Please adjust len_correration and/or new_time_origin."
            )
            return
        if 2 * corr_steps - 1 >= self._num_frame:
            print(
                "VACF correlation time is too long for the number of frames. "
                "Please reduce len_correration or use a longer trajectory."
            )
            return

        if self._is_diffusion_bin:
            print("Binning and VACF-approaches cannot be run in parallel.")
            return
        if self._is_diffusion_mc:
            print("MC and VACF-approaches cannot be run in parallel.")
            return
        if self._is_numpy:
            print("VACF and numpy-approaches cannot be run in parallel.")
            return
        if direction not in [0, 1, 2] and not (direction == "radial" and self._pore):
            print(
                "Wrong directional input. Options: 0 (x), 1 (y), 2 (z), "
                "or 'radial' (pore systems only)."
            )
            return
        if self._traj.split(".")[-1] != "trr":
            print("VACF requires a .trr trajectory file with velocities.")
            return

        if direction == "radial":
            axes, shapes = [], []
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    if self._pore_props[pore_id]["type"] == "CYLINDER":
                        axes.append(tuple(self._pore_props[pore_id]["focal"][:2]))
                        shapes.append(self._pore_props[pore_id]["type"])
                    else:
                        print(
                            "Radial VACF only available for cylindrical pore systems."
                        )
                        return
            if not all(a == axes[0] for a in axes):
                print(
                    "Radial VACF requires all pores to share the same cylindrical axis."
                )
                return
            if all(s == shapes[0] for s in shapes) and shapes[0] == "CYLINDER":
                direction = "radial_cylindrical"
            else:
                print("Radial VACF requires all pores to be the same shape.")
                return

        self._is_diffusion_vacf = True

        if direction == "radial_cylindrical":
            bins = self._bin_pore(bin_num)["bins"]
        else:
            bins = self._bin_mc(bin_num, direction)["bins"]

        self._diff_vacf_inp = {
            "output": link_out,
            "bins": bins,
            "len_correration": len_correration,
            "new_time_origin": new_time_origin,
            "sample_step": sample_step,
            "len_frame": len_frame,
            "bin_num": bin_num,
            "direction": direction,
            "corr_steps": corr_steps,
            "new_time_origin_steps": new_time_origin_steps,
            "num_res": self.num_res,
            "sample_each_residue": sample_each_residue,
        }

    def _diffusion_vacf_data(self) -> dict:
        """Create the VACF diffusion data structure."""
        bin_num = self._diff_vacf_inp["bin_num"]
        corr_steps = self._diff_vacf_inp["corr_steps"]
        per_res = self._diff_vacf_inp["sample_each_residue"]
        n_res = self.num_res if per_res else 1

        data = {}
        if self._pore:
            for pore_id in self._pore:
                if pore_id[:5] == "shape":
                    data[pore_id] = {
                        "vacf_data": np.zeros((bin_num, corr_steps, n_res, 3), float),
                        "density": np.zeros((bin_num, n_res), int),
                    }
        else:
            data = {
                "vacf_data": np.zeros((bin_num, corr_steps, n_res, 3), float),
                "density": np.zeros((bin_num, n_res), int),
            }
        return data

    def _diffusion_vacf(
        self,
        data: dict,
        frame_id: int,
        pos_list: np.ndarray,
        pos_pointer: int,
        vel_list: np.ndarray,
        vel_pointer: int,
    ):
        """Sample the local VACF for all bins at one time origin.

        Parameters
        ----------
        data : dictionary
            VACF data structure
        frame_id : integer
            Current frame index
        pos_list : numpy.ndarray
            Circular buffer of CoM positions, shape (corr_steps, n_res, 3)
        pos_pointer : integer
            Current write position in pos_list
        vel_list : numpy.ndarray
            Circular buffer of CoM velocities, shape (2*corr_steps-1, n_res, 3)
        vel_pointer : integer
            Current write position in vel_list
        """
        bins = self._diff_vacf_inp["bins"]
        direction = self._diff_vacf_inp["direction"]
        corr_steps = self._diff_vacf_inp["corr_steps"]
        new_time_origin_steps = self._diff_vacf_inp["new_time_origin_steps"]
        per_res = self._diff_vacf_inp["sample_each_residue"]

        if direction == "radial_cylindrical":
            direction = 0

        if (frame_id - corr_steps) % new_time_origin_steps != 0:
            return

        pos = pos_list[pos_pointer, :, direction]
        vel_pointer = (vel_pointer + corr_steps - 1) % vel_list.shape[0]
        forward_idx = (np.arange(corr_steps) + vel_pointer) % vel_list.shape[0]
        backward_idx = (
            np.arange(corr_steps)[::-1] + vel_pointer + corr_steps
        ) % vel_list.shape[0]

        for pore_id in self._pore.keys() if self._pore else [None]:
            if pore_id is not None:
                if pore_id[:5] != "shape":
                    continue
                data_p = data[pore_id]
                bin_p = bins[pore_id]
                res = self._pore_props["box"]["res"]
                z_min = (
                    res
                    + self._pore_props[pore_id]["focal"][2]
                    - self._pore_props[pore_id]["length"] / 2
                    + self._entry
                )
                z_max = (
                    res
                    + self._pore_props[pore_id]["focal"][2]
                    + self._pore_props[pore_id]["length"] / 2
                    - self._entry
                )
                pore_mask = (z_min < pos_list[pos_pointer, :, 2]) & (
                    pos_list[pos_pointer, :, 2] < z_max
                )
                in_wall_mask = pos > self._pore_props[pore_id]["diam"] * 1.01 / 2
            else:
                data_p = data
                bin_p = bins
                pore_mask = np.ones(pos.shape, dtype=bool)
                in_wall_mask = np.zeros(pos.shape, dtype=bool)

            com_bins = np.digitize(pos, bin_p) - 1
            n_bins = len(bin_p) - 1
            for bin_id in range(n_bins):
                mask = (com_bins == bin_id) & pore_mask & ~in_wall_mask
                if not np.any(mask):
                    continue
                vel0 = vel_list[vel_pointer, mask, :]
                fwd = vel_list[forward_idx][:, mask, :]
                bwd = vel_list[backward_idx][:, mask, :]
                vacf = (vel0[np.newaxis] * fwd + vel0[np.newaxis] * bwd) / 2

                if per_res:
                    idx = np.where(mask)[0]
                    data_p["vacf_data"][bin_id, :, idx, :] += vacf.transpose(0, 1, 2)[
                        :, :, :
                    ]
                    data_p["density"][bin_id, idx] += 1
                else:
                    data_p["vacf_data"][bin_id, :, 0, :] += np.sum(vacf, axis=1)
                    data_p["density"][bin_id, 0] += np.sum(mask)

    ##################
    # Numpy Sampling #
    ##################
    def init_numpy_file(self, link_out, positions=True, velocities=True):
        """Enable numpy position/velocity sampling.

        Parameters
        ----------
        link_out : string
            Output ``.npz`` file path
        positions : bool, optional
            Sample CoM positions (default True)
        velocities : bool, optional
            Sample CoM velocities (default True); requires a ``.trr`` trajectory
        """
        if self._is_diffusion_bin or self._is_diffusion_mc or self._is_diffusion_vacf:
            print("Numpy sampling cannot run in parallel with diffusion sampling.")
            return
        if velocities and self._traj.split(".")[-1] != "trr":
            print("Velocity sampling requires a .trr trajectory file.")
            return
        self._is_numpy = True
        self._numpy_inp = {
            "output": link_out,
            "positions": positions,
            "velocities": velocities,
        }

    def _numpy_data(self):
        """Create numpy sampling data structure."""
        data = {}
        if self._numpy_inp["positions"]:
            data["positions"] = np.zeros((self._num_frame, self.num_res, 3), float)
        if self._numpy_inp["velocities"]:
            data["velocities"] = np.zeros((self._num_frame, self.num_res, 3), float)
        return data

    def _numpy(self, data, positions, velocities, frame_id):
        """Record CoM positions and velocities for one frame.

        Parameters
        ----------
        data : dictionary
            Numpy data structure
        positions : array-like
            Raw atom positions for this frame
        velocities : array-like
            Raw atom velocities for this frame
        frame_id : integer
            Current frame index
        """
        if self._numpy_inp["positions"]:
            pos = (
                np.asarray(positions).reshape((self.num_res, self._atoms_per_mol, 3))[
                    :, self._atoms, :
                ]
                / 10
            )
            data["positions"][frame_id] = (
                np.sum(pos * self._masses_arr[np.newaxis, :, np.newaxis], axis=1)
                / self._sum_masses
            )
        if self._numpy_inp["velocities"]:
            vel = (
                np.asarray(velocities).reshape((self.num_res, self._atoms_per_mol, 3))[
                    :, self._atoms, :
                ]
                * 100
            )
            data["velocities"][frame_id] = (
                np.sum(vel * self._masses_arr[np.newaxis, :, np.newaxis], axis=1)
                / self._sum_masses
            )

    ############
    # Sampling #
    ############
    def sample(
        self, shift=[0, 0, 0], n_proc=0, is_pbc=True, is_broken=False, is_parallel=True
    ):
        """Run all enabled sampling routines.

        Parameters
        ----------
        shift : list, optional
            Translation vector in nm
        n_proc : integer, optional
            Number of CPU cores (0 = all available)
        is_pbc : bool, optional
            True to apply periodic boundary conditions
        is_broken : bool, optional
            True to flag broken molecules
        is_parallel : bool, optional
            True to use multiprocessing
        """
        if len(shift) != 3:
            print("Sample - Wrong shift dimension.")
            return

        n_proc = n_proc if n_proc and n_proc <= mp.cpu_count() else mp.cpu_count()

        if is_parallel and self._is_angle:
            print("Currently the angle routine cannot be parallelized...")
            return

        if is_parallel:
            frame_num = math.floor(self._num_frame / n_proc)
            frame_start = [frame_num * i for i in range(n_proc)]
            frame_end = [
                frame_num * (i + 1) if i < n_proc - 1 else self._num_frame
                for i in range(n_proc)
            ]

            if self._is_diffusion_bin:
                win = self._diff_bin_inp["len_window"] * self._diff_bin_inp["len_step"]
                frame_start = [
                    x - win + 1 if i > 0 else x for i, x in enumerate(frame_start)
                ]

            if self._is_diffusion_mc:
                max_step = max(self._diff_mc_inp["len_step"])
                frame_end = [x + max_step for x in frame_end]
                for i in range(len(frame_end)):
                    if frame_end[i] >= self._num_frame:
                        frame_end[i] = frame_end[-1] - max_step

            if self._is_diffusion_vacf:
                cs = self._diff_vacf_inp["corr_steps"]
                nto = self._diff_vacf_inp["new_time_origin_steps"]
                frame_start = [max(0, s - s % nto - cs + 1) for s in frame_start]
                frame_end = [min(self._num_frame, e - e % nto + cs) for e in frame_end]

            frame_np = [
                list(range(frame_start[i], frame_end[i])) for i in range(n_proc)
            ]
            _ctx = mp.get_context("fork") if sys.platform != "win32" else mp
            pool = _ctx.Pool(processes=n_proc)
            results = [
                pool.apply_async(
                    self._sample_helper, args=(fl, shift, is_pbc, is_broken)
                )
                for fl in frame_np
            ]
            pool.close()
            pool.join()
            output = [x.get() for x in results]
            del results
        else:
            output = [
                self._sample_helper(
                    list(range(self._num_frame)), shift, is_pbc, is_broken
                )
            ]

        system = (
            {"sys": "pore", "props": self._pore_props}
            if self._pore
            else {"sys": "box", "props": {"length": self._box}}
        )
        inp = {
            "num_frame": self._num_frame,
            "mass": self._mol.get_mass(),
            "entry": self._entry,
        }

        if self._is_density:
            inp_dens = {
                **inp,
                **{k: v for k, v in self._dens_inp.items() if k != "output"},
            }
            data_dens = output[0]["density"]
            for out in output[1:]:
                if self._pore:
                    for pore_id in data_dens:
                        if pore_id[:5] == "shape":
                            data_dens[pore_id]["in"] = [
                                x + y
                                for x, y in zip(
                                    data_dens[pore_id]["in"],
                                    out["density"][pore_id]["in"],
                                )
                            ]
                data_dens["ex"] = [
                    x + y for x, y in zip(data_dens["ex"], out["density"]["ex"])
                ]
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_dens,
                    "data": data_dens,
                    "type": "dens_bin",
                },
                self._dens_inp["output"],
            )

        if self._is_gyration:
            inp_gyr = {
                **inp,
                **{k: v for k, v in self._gyr_inp.items() if k != "output"},
            }
            data_gyr = output[0]["gyration"]
            for out in output[1:]:
                if self._pore:
                    for pore_id in data_gyr:
                        if pore_id[:5] == "shape":
                            data_gyr[pore_id]["in"] = [
                                x + y
                                for x, y in zip(
                                    data_gyr[pore_id]["in"],
                                    out["gyration"][pore_id]["in"],
                                )
                            ]
                data_gyr["ex"] = [
                    x + y for x, y in zip(data_gyr["ex"], out["gyration"]["ex"])
                ]
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_gyr,
                    "data": data_gyr,
                    "type": "gyr_bin",
                },
                self._gyr_inp["output"],
            )

        if self._is_angle:
            inp_angle = {
                **inp,
                **{k: v for k, v in self._angle_inp.items() if k != "output"},
            }
            data_angle = output[0]["angle"]
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_angle,
                    "data": data_angle,
                    "type": "angle_bin",
                },
                self._angle_inp["output"],
            )

        if self._is_diffusion_bin:
            inp_diff = {
                **inp,
                **{k: v for k, v in self._diff_bin_inp.items() if k != "output"},
            }
            data_diff = output[0]["diffusion_bin"]
            for out in output[1:]:
                for pore_id in data_diff:
                    if pore_id[:5] == "shape":
                        for key in ("z", "r", "n", "z_tot", "r_tot", "n_tot"):
                            a = np.array(data_diff[pore_id][key])
                            b = np.array(out["diffusion_bin"][pore_id][key])
                            data_diff[pore_id][key] = (a + b).tolist()
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_diff,
                    "data": data_diff,
                    "type": "diff_bin",
                },
                self._diff_bin_inp["output"],
            )

        if self._is_diffusion_mc:
            inp_diff = {
                **inp,
                **{k: v for k, v in self._diff_mc_inp.items() if k != "output"},
            }
            data_diff = output[0]["diffusion_mc"]
            for step in self._diff_mc_inp["len_step"]:
                for out in output[1:]:
                    data_diff[step] += out["diffusion_mc"][step]
            for step in self._diff_mc_inp["len_step"]:
                data_diff[step] = data_diff[step][1:-1, 1:-1]
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_diff,
                    "data": data_diff,
                    "type": "diff_mc",
                },
                self._diff_mc_inp["output"],
            )

        if self._is_diffusion_vacf:
            inp_diff = {
                **inp,
                **{k: v for k, v in self._diff_vacf_inp.items() if k != "output"},
            }
            data_diff = output[0]["diffusion_vacf"]
            for out in output[1:]:
                for pore_id in self._pore.keys() if self._pore else [None]:
                    if pore_id is not None:
                        if pore_id[:5] != "shape":
                            continue
                        data_diff[pore_id]["density"] += out["diffusion_vacf"][pore_id][
                            "density"
                        ]
                        data_diff[pore_id]["vacf_data"] += out["diffusion_vacf"][
                            pore_id
                        ]["vacf_data"]
                    else:
                        data_diff["density"] += out["diffusion_vacf"]["density"]
                        data_diff["vacf_data"] += out["diffusion_vacf"]["vacf_data"]
            utils.save(
                {
                    system["sys"]: system["props"],
                    "inp": inp_diff,
                    "data": data_diff,
                    "type": "diff_vacf",
                },
                self._diff_vacf_inp["output"],
            )

        if self._is_numpy:
            data_np = output[0]["numpy"].copy()
            for out in output[1:]:
                if self._numpy_inp["positions"]:
                    data_np["positions"] += out["numpy"]["positions"]
                if self._numpy_inp["velocities"]:
                    data_np["velocities"] += out["numpy"]["velocities"]
            np.savez(self._numpy_inp["output"], **data_np)

    def _sample_helper(self, frame_list, shift, is_pbc, is_broken):
        """Worker function: sample all enabled routines for a list of frames.

        Parameters
        ----------
        frame_list : list
            Frame ids to process
        shift : list
            Translation vector in nm
        is_pbc : bool
            Apply periodic boundary conditions
        is_broken : bool
            Check for broken molecules

        Returns
        -------
        output : dictionary
            Accumulated sampling data
        """
        box = (
            np.array(self._pore_props["box"]["dimensions"])
            if self._pore
            else np.array(self._box)
        )
        res = self._pore_props["box"]["res"] if self._pore else 0
        shift_arr = np.array(shift, dtype=float)

        # Initialize sliding windows for diffusion
        if self._is_diffusion_bin:
            com_list = {pid: [] for pid in self._pore if pid[:5] == "shape"}
            idx_list = {pid: [] for pid in self._pore if pid[:5] == "shape"}
        elif self._is_diffusion_mc:
            com_list = []
            idx_list = []

        # VACF circular buffers
        if self._is_diffusion_vacf:
            cs = self._diff_vacf_inp["corr_steps"]
            _pos_list = np.zeros((cs, self.num_res, 3), float)
            _pos_ptr = 0
            _vel_list = np.zeros((2 * cs - 1, self.num_res, 3), float)
            _vel_ptr = 0
            _vacf_filled = False

        # Initialize per-worker output structures
        output = {}
        if self._is_density:
            output["density"] = self._density_data()
        if self._is_gyration:
            output["gyration"] = self._gyration_data()
        if self._is_angle:
            output["angle"] = self._angle_data()
        if self._is_diffusion_bin:
            output["diffusion_bin"] = self._diffusion_bin_data()
        if self._is_diffusion_mc:
            output["diffusion_mc"] = self._diffusion_mc_data()
        if self._is_diffusion_vacf:
            output["diffusion_vacf"] = self._diffusion_vacf_data()
        if self._is_numpy:
            output["numpy"] = self._numpy_data()

        # Sliding window fill length
        if self._is_diffusion_bin:
            len_fill = self._diff_bin_inp["len_window"] * self._diff_bin_inp["len_step"]
        elif self._is_diffusion_mc:
            len_fill = self._diff_mc_inp["len_step"][-1] + 1
        else:
            len_fill = 1

        traj = cf.Trajectory(self._traj)
        frame_fmt = "%{}i".format(len(str(self._num_frame)))

        for frame_id in frame_list:
            frame = traj.read_step(frame_id)
            positions = frame.positions  # shape (N, 3) in Angstroms

            # ── VACF / numpy batch processing (all residues at once) ──────────
            if self._is_diffusion_vacf or self._is_numpy:
                all_pos = (
                    np.asarray(positions).reshape(self.num_res, self._atoms_per_mol, 3)[
                        :, self._atoms, :
                    ]
                    / 10
                    + shift_arr
                )  # (n_res, n_sel, 3) in nm
                all_com = (
                    np.sum(
                        all_pos * self._masses_arr[np.newaxis, :, np.newaxis], axis=1
                    )
                    / self._sum_masses
                )  # (n_res, 3)
                if is_pbc:
                    all_com = all_com - np.floor(all_com / box) * box

                if self._is_diffusion_vacf:
                    all_vel = (
                        np.asarray(frame.velocities).reshape(
                            self.num_res, self._atoms_per_mol, 3
                        )[:, self._atoms, :]
                        * 100
                    )  # Å/ps → m/s
                    all_vel_com = (
                        np.sum(
                            all_vel * self._masses_arr[np.newaxis, :, np.newaxis],
                            axis=1,
                        )
                        / self._sum_masses
                    )

                    _com_vacf = all_com.copy()
                    _vel_vacf = all_vel_com.copy()
                    if self._diff_vacf_inp["direction"] == "radial_cylindrical":
                        focal = np.array(self._pore_props["shape_00"]["focal"][:2])
                        dx = all_com[:, 0] - focal[0]
                        dy = all_com[:, 1] - focal[1]
                        r = np.sqrt(dx**2 + dy**2)
                        r_safe = np.where(r == 0, 1e-8, r)
                        _com_vacf[:, 0] = r
                        _com_vacf[:, 1] = np.arctan2(dy, dx)
                        vr = (all_vel_com[:, 0] * dx + all_vel_com[:, 1] * dy) / r_safe
                        vt = (all_vel_com[:, 1] * dx - all_vel_com[:, 0] * dy) / r_safe
                        _vel_vacf[:, 0] = vr
                        _vel_vacf[:, 1] = vt

                    _pos_list[_pos_ptr] = _com_vacf
                    _pos_ptr = (_pos_ptr + 1) % cs
                    _vel_list[_vel_ptr] = _vel_vacf
                    _vel_ptr += 1
                    if _vel_ptr >= 2 * cs - 1:
                        _vel_ptr = 0
                        _vacf_filled = True

                    if _vacf_filled:
                        self._diffusion_vacf(
                            output["diffusion_vacf"],
                            frame_id,
                            _pos_list,
                            _pos_ptr,
                            _vel_list,
                            _vel_ptr,
                        )

                if self._is_numpy:
                    self._numpy(output["numpy"], positions, frame.velocities, frame_id)

            # Manage sliding window lists
            if self._is_diffusion_bin:
                for pid in com_list:
                    if len(com_list[pid]) >= len_fill:
                        idx_list[pid].pop(0)
                        com_list[pid].pop(0)
                    idx_list[pid].append({})
                    com_list[pid].append({})
            elif self._is_diffusion_mc:
                if len(com_list) >= len_fill:
                    idx_list.pop(0)
                    com_list.pop(0)
                idx_list.append({})
                com_list.append({})

            for res_id, atom_indices in self._res_list.items():
                # Extract and convert positions: Å → nm, apply shift
                pos = positions[atom_indices] / 10.0 + shift_arr  # (n_atoms, 3)

                # Centre of mass without PBC
                com_no_pbc = (self._masses_arr @ pos) / self._sum_masses  # (3,)

                if is_broken and np.any(np.abs(com_no_pbc - pos[0]) > box / 3):
                    print(f"Sample - Broken molecule found - ResID: {res_id:5d}")

                com = com_no_pbc % box if is_pbc else com_no_pbc

                # Distance from pore axis
                if self._pore:
                    dist = {}
                    for pore_id, pp in self._pore_props.items():
                        if pore_id[:5] == "shape":
                            focal = pp["focal"]
                            if pp["type"] in ("CYLINDER", "CONE"):
                                dist[pore_id] = float(
                                    np.linalg.norm(
                                        com[:2] - np.array([focal[0], focal[1]])
                                    )
                                )
                            elif pp["type"] == "SLIT":
                                dist[pore_id] = abs(focal[1] - float(com[1]))
                else:
                    dist = 0

                # Region classification
                region = ""
                pore_in = 1
                if (
                    self._pore
                    and res + self._entry < com[2] < box[2] - res - self._entry
                ):
                    region = "in"
                    for pore_id, pp in self._pore_props.items():
                        if pore_id[:5] == "shape":
                            z_min = (
                                res + pp["focal"][2] - pp["length"] / 2 + self._entry
                            )
                            z_max = (
                                res + pp["focal"][2] + pp["length"] / 2 - self._entry
                            )
                            if (
                                z_min < com[2] < z_max
                                and dist[pore_id] < pp["diam"] * 1.01 / 2
                            ):
                                pore_in = pore_id
                elif not self._pore or com[2] < res or com[2] > box[2] - res:
                    region = "ex"
                    pore_in = 0

                # Determine whether to sample (skip window-filling frames for non-first workers)
                if self._is_diffusion_bin:
                    is_sample = any(
                        len(com_list[pid]) == len_fill or frame_id <= len_fill
                        for pid in com_list
                    )
                else:
                    is_sample = True

                if is_sample:
                    if self._is_density and pore_in != 1:
                        self._density(output["density"], region, dist, com, pore_in)
                    if self._is_gyration and pore_in != 1:
                        self._gyration(
                            output["gyration"], region, dist, com_no_pbc, pos, pore_in
                        )
                    if self._is_angle and pore_in != 1:
                        self._angle(output["angle"], region, dist, com, pos, pore_in)

                if self._is_diffusion_bin and pore_in != 1:
                    self._diffusion_bin(
                        output["diffusion_bin"],
                        region,
                        pore_in,
                        dist,
                        com_list,
                        idx_list,
                        res_id,
                        com,
                    )
                if self._is_diffusion_mc:
                    self._diffusion_mc(
                        output["diffusion_mc"],
                        idx_list,
                        com,
                        res_id,
                        frame_list,
                        frame_id,
                    )

            if (
                (frame_id + 1) % 10 == 0
                or frame_id == 0
                or frame_id == self._num_frame - 1
            ):
                sys.stdout.write(
                    "Finished frame "
                    + frame_fmt % (frame_id + 1)
                    + "/"
                    + frame_fmt % self._num_frame
                    + "...\r"
                )
                sys.stdout.flush()
        print()
        return output
