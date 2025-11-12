.. container:: cell markdown

   .. rubric:: Diffusion Analysis with the VACF
      :name: diffusion-analysis-with-the-vacf

   Example code to analyse a local diffusion profile over the simulation
   box. The method can run for a specific pore system (using a pore.yml
   file from PoreMS) or an other simulation box.

   .. rubric:: Import packages
      :name: import-packages

.. container:: cell code

   .. code:: python

      import poreana as pa
      import porems as pms
      import matplotlib.pyplot as plt

.. container:: cell code

   .. code:: python

      # Plot settings seaborn
      import seaborn as sns
      sns.set(font_scale=5) 
      sns.set_context("poster")
      plt.rc('font', size = 25) # steuert die Standardtextgröße
      plt.rc('axes', titlesize = 25) # Schriftgröße des Titels
      plt.rc('axes', labelsize = 25) # Schriftgröße der x- und y-Beschriftungen
      plt.rc('xtick', labelsize = 25) #Schriftgröße der x-Tick-Labels
      plt.rc('ytick', labelsize = 25) #Schriftgröße der y-Tick-Labels
      plt.rc('legend', fontsize = 25) #Schriftgröße der Legende
      sns.set_style("white", {"xtick.bottom": True,'ytick.left': True})

.. container:: cell markdown

   .. rubric:: Create the molecules
      :name: create-the-molecules

   Load the molecule to be analyzed (here: water, using the TIP4P/2005
   force field). \\ Because this force field defines a specific dipole
   moment, we must explicitly set the atomic masses. \\ For most other
   molecules, this step is not required.

.. container:: cell code

   .. code:: python

      # Load molecule
      mol = pms.Molecule(inp="tip4p2005.gro")

      # Set mass
      mol.set_masses([15.999,1.00784,1.00784,0])

.. container:: cell markdown

   .. rubric:: VACF Inputs
      :name: vacf-inputs

   Specifiy VACF inputs for the sampling

.. container:: cell code

   .. code:: python

      # Inputs for sampling
      len_corr = 15e-12           # Correlation time of the molecule in seconds
      new_time_o = 1e-13          # New time origin in seconds
      sample_st = 20              # sample rate of frames
      len_fr = 1e-15              # time between two frames in seconds

.. container:: cell markdown

   .. rubric:: Sample the VACF for the trajectory
      :name: sample-the-vacf-for-the-trajectory

   The local VACF diffusion analysis need the sampled object file using
   the vacf diffusion routine. In this example the density is sampled
   along the z-axis (direction = 2) of the simulation box.

.. container:: cell code

   .. code:: python

      sample = pa.Sample("pore.yml", "traj_SOL.trr", mol, masses = [15.999,1.00784,1.00784,0])
      sample.init_diffusion_vacf("diff_vacf2.obj", len_correration=len_corr, new_time_origin=new_time_o, sample_step=sample_st, len_frame=len_fr, bin_num=32, direction=2)
      sample.sample(is_parallel=False)

   .. container:: output stream stderr

      ::

         /home/marc/Nextcloud/codes/PoreAna/poreana/sample.py:92: ChemfilesWarning: could not open the file at 'traj_SOL.trr'
           traj = cf.Trajectory(self._traj)

   .. container:: output error

      ::

         ---------------------------------------------------------------------------
         ValueError                                Traceback (most recent call last)
         File ~/anaconda3/lib/python3.11/site-packages/chemfiles/utils.py:102, in _check_handle(handle)
             101 try:
         --> 102     handle.contents
             103 except ValueError:

         ValueError: NULL pointer access

         During handling of the above exception, another exception occurred:

         ChemfilesError                            Traceback (most recent call last)
         Cell In[4], line 1
         ----> 1 sample = pa.Sample("pore.yml", "traj_SOL.trr", mol, masses = [15.999,1.00784,1.00784,0])
               2 sample.init_diffusion_vacf("diff_vacf2.obj", len_correration=len_corr, new_time_origin=new_time_o, sample_step=sample_st, len_frame=len_fr, bin_num=32, direction=2)
               3 sample.sample(is_parallel=False)

         File ~/Nextcloud/codes/PoreAna/poreana/sample.py:92, in Sample.__init__(self, system, link_traj, mol, atoms, masses, entry, frame_end)
              89     return
              91 # Get number of frames
         ---> 92 traj = cf.Trajectory(self._traj)
              93 if frame_end == -1:
              94     self._num_frame = traj.nsteps

         File ~/anaconda3/lib/python3.11/site-packages/chemfiles/trajectory.py:155, in Trajectory.__init__(self, path, mode, format)
             153 self.__mode = mode
             154 self.__format = format
         --> 155 super(Trajectory, self).__init__(ptr)

         File ~/anaconda3/lib/python3.11/site-packages/chemfiles/trajectory.py:22, in BaseTrajectory.__init__(self, ptr)
              20 def __init__(self, ptr):
              21     self.__closed = False
         ---> 22     super(BaseTrajectory, self).__init__(ptr, is_const=False)

         File ~/anaconda3/lib/python3.11/site-packages/chemfiles/utils.py:35, in CxxPointer.__init__(self, ptr, is_const, origin)
              33 self.__origin = origin
              34 self.__frozen = True
         ---> 35 _check_handle(ptr)

         File ~/anaconda3/lib/python3.11/site-packages/chemfiles/utils.py:104, in _check_handle(handle)
             102     handle.contents
             103 except ValueError:
         --> 104     raise ChemfilesError(_last_error())

         ChemfilesError: could not open the file at 'traj_SOL.trr'

.. container:: cell markdown

   The VACF is used to sample the diffusion inside a pore in radial,
   axial and tangial direction (direction = "radial").

.. container:: cell code

   .. code:: python

      sample = pa.Sample("pore.yml", "traj_SOL.trr", mol, masses = [15.999,1.00784,1.00784,0])
      sample.init_diffusion_vacf("diff_vacf_radial.obj", len_correration=len_corr, new_time_origin=new_time_o, sample_step=sample_st, len_frame=len_fr, bin_num=32, direction="radial")
      sample.sample(is_parallel=False)

.. container:: cell markdown

   .. rubric:: Display the integrated VACF for every bin
      :name: display-the-integrated-vacf-for-every-bin

   With the sampling obj-file the integrated velocity correlation for
   every bin in the simulation box can be displayed.

.. container:: cell code

   .. code:: python

      fig, ax = plt.subplots(figsize=(30,9))
      ax1 = plt.subplot(1,2,1)
      plt.title("Box sampling")
      pa.diffusion.plot_correlation_per_bin("diff_vacf.obj", plot_axis = ax1, direction= "z")
      ax2 = plt.subplot(1,2,2)
      plt.title("Pore sampling")
      pa.diffusion.plot_correlation_per_bin("diff_vacf_radial.obj", plot_axis = ax2, pore_id = "shape_00", direction = "a")
      plt.xlim([0,2])
      plt.tight_layout()
      plt.savefig("integration_vacf.svg")

   .. container:: output stream stderr

      ::

         /tmp/ipykernel_1997484/3463115208.py:2: MatplotlibDeprecationWarning: Auto-removal of overlapping axes is deprecated since 3.6 and will be removed two minor releases later; explicitly call ax.remove() as needed.
           ax1 = plt.subplot(1,2,1)

   .. container:: output stream stdout

      ::

         Sampled 673_047_282 data points (including time reversal) for VACF calculation.
         Sampled 145_760_688 data points (including time reversal) for VACF calculation.

   .. container:: output display_data

      |image1|

.. container:: cell markdown

   .. rubric:: Display the diffusion profile over the box length
      :name: display-the-diffusion-profile-over-the-box-length

   The following function can use to show the diffusion profile over the
   simulation box

.. container:: cell code

   .. code:: python

      fig, ax = plt.subplots(figsize=(30,9))
      ax1 = plt.subplot(1,2,1)
      plt.title("Box sampling")
      diffusion_profile, diffusion_mean = pa.diffusion.diffusion_per_bin("diff_vacf.obj", plot_selection="xyz", plot_axis= ax1, combine_bins=4, mean_over_time=10e-12)
      plt.legend()
      ax2 = plt.subplot(1,2,2)
      plt.title("Pore sampling")
      diffusion_profile, diffusion_mean = pa.diffusion.diffusion_per_bin("diff_vacf_radial.obj", plot_selection="rat", pore_id = "shape_00", plot_axis= ax2, combine_bins=6, mean_over_time=10e-12)
      plt.legend()
      plt.savefig("diff_profile_vacf.svg")

   .. container:: output stream stderr

      ::

         /tmp/ipykernel_1997484/3804608588.py:2: MatplotlibDeprecationWarning: Auto-removal of overlapping axes is deprecated since 3.6 and will be removed two minor releases later; explicitly call ax.remove() as needed.
           ax1 = plt.subplot(1,2,1)

   .. container:: output stream stdout

      ::

         Sampled 673_047_282 data points (including time reversal) for VACF calculation.
         Mean over last 500 steps.
         Sampled 145_760_688 data points (including time reversal) for VACF calculation.
         Mean over last 500 steps.

   .. container:: output display_data

      |image2|

.. container:: cell markdown

   .. rubric:: Calculating diffusion coefficients
      :name: calculating-diffusion-coefficients

.. container:: cell code

   .. code:: python

      # Using the diffusion profile in z-directions
      diffusion_profile, diffusion_mean_res = pa.diffusion.diffusion_per_bin("diff_vacf.obj", section = [0,5], combine_bins=8, mean_over_time=10e-12)
      diffusion_profile, diffusion_mean_pore = pa.diffusion.diffusion_per_bin("diff_vacf.obj", section = [5,15], combine_bins=4, mean_over_time=10e-12)
      print("\nReservoir diffusion : ", diffusion_mean_res[2], "m²/s")
      print("Pore diffusion : ", diffusion_mean_pore[2], "m²/s")
      print("Reservoir/Pore ratio : ", diffusion_mean_res[2]/diffusion_mean_pore[2])

      # Using the diffusion profile in the pore
      diffusion_profile, diffusion_rad_pore = pa.diffusion.diffusion_per_bin("diff_vacf_radial.obj",pore_id = "shape_00", combine_bins=4, mean_over_time=10e-12)
      print("Pore diffusion (axial) : ", diffusion_rad_pore[2], "m²/s")
      print("Pore diffusion (radial) : ", diffusion_rad_pore[0], "m²/s")

   .. container:: output stream stdout

      ::

         Sampled 673_047_282 data points (including time reversal) for VACF calculation.
         Mean over last 500 steps.
         Sampled 673_047_282 data points (including time reversal) for VACF calculation.
         Mean over last 500 steps.

         Reservoir diffusion :  1.9675621308108457 m²/s
         Pore diffusion :  1.457195798962518 m²/s
         Reservoir/Pore ratio :  1.3502386791203311
         Sampled 145_760_688 data points (including time reversal) for VACF calculation.
         Mean over last 500 steps.
         Pore diffusion (axial) :  1.5062768636353063 m²/s
         Pore diffusion (radial) :  1.2449438793219074 m²/s

.. |image1| image:: af023bfaf77e137d892e65b044c544873d618616.png
.. |image2| image:: 6c71d6954c099166a33abf57064b5e27e8652acd.png
