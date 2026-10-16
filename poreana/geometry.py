################################################################################
# Geometry                                                                     #
#                                                                              #
"""Here basic geometric functions are noted."""
################################################################################

import numpy as np


def dot_product(vec_a, vec_b):
    r"""Calculate the dot product of two vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{dot}(\boldsymbol{a},\boldsymbol{b})=
        \begin{pmatrix}a_1\\\vdots\\a_n\end{pmatrix}\cdot
        \begin{pmatrix}b_1\\\vdots\\b_n\end{pmatrix}=
        a_1\cdot b_1+a_2\cdot b_2+\dots+a_n\cdot b_n.

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`

    Returns
    -------
    dot : float
        Dot product value
    """
    return float(np.dot(vec_a, vec_b))


def length(vec):
    r"""Calculate the length of a vector
    :math:`\boldsymbol{a}\in\mathbb{R}^n`

    .. math::

        \text{length}(\boldsymbol{a})=|\boldsymbol{a}|
        =\sqrt{\boldsymbol{a}\cdot\boldsymbol{a}}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    length : float
        Vector length
    """
    return float(np.linalg.norm(vec))


def vector(pos_a, pos_b):
    r"""Calculate the vector between to two positions
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{vec}(\boldsymbol{a},\boldsymbol{b})
        =\begin{pmatrix}b_1-a_1\\\vdots\\b_n-a_n\end{pmatrix}

    Parameters
    ----------
    pos_a : list
        First position :math:`\boldsymbol{a}`
    pos_b : list
        Second position :math:`\boldsymbol{b}`

    Returns
    -------
    vector : list
        Bond vector
    """
    if len(pos_a) != len(pos_b):
        print("Vector: Wrong dimensions...")
        return

    return list(np.asarray(pos_b) - np.asarray(pos_a))


def unit(vec):
    r"""Transform a vector :math:`\boldsymbol{a}\in\mathbb{R}^n` into a
    unit vector

    .. math::

        \text{unit}(\boldsymbol{a})
        =\frac{\boldsymbol{a}}{|\boldsymbol{a}|}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    vec : list
        Unit vector
    """
    v = np.asarray(vec, dtype=float)
    n = np.linalg.norm(v)
    return list(v / n if n != 0 else v)


def cross_product(vec_a, vec_b):
    r"""Calculate the cross product of two three-dimensional vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^3`

    .. math::

        \text{cross}(\boldsymbol{a},\boldsymbol{b})=\begin{pmatrix}
        a_2\cdot b_3-a_3\cdot b_2\\
        a_3\cdot b_1-a_1\cdot b_4\\
        a_1\cdot b_2-a_2\cdot b_1
        \end{pmatrix}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`

    Returns
    -------
    vec : list
        Cross product vector
    """
    return list(np.cross(vec_a, vec_b))


def angle(vec_a, vec_b, is_deg=True):
    r"""Calculate the angle between two vectors
    :math:`\boldsymbol{a},\boldsymbol{b}\in\mathbb{R}^n`

    .. math::

        \text{angle}=\cos^{-1}\frac{\boldsymbol{a}\cdot\boldsymbol{b}}
        {|\boldsymbol{a}||\boldsymbol{a}|}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\boldsymbol{b}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Angle
    """
    a = np.asarray(vec_a, dtype=float)
    b = np.asarray(vec_b, dtype=float)
    cos_val = np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1.0, 1.0)
    rad = np.arccos(cos_val)
    return float(np.degrees(rad) if is_deg else rad)
