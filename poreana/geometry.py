################################################################################
# Geometry                                                                     #
#                                                                              #
"""Here basic geometric functions are noted."""
################################################################################


import math
from operator import mul


def dot_product(vec_a, vec_b):
    """Calculate the dot product of two vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{dot}(\\boldsymbol{a},\\boldsymbol{b})=
        \\begin{pmatrix}a_1\\\\\\vdots\\\\a_n\\end{pmatrix}\\cdot
        \\begin{pmatrix}b_1\\\\\\vdots\\\\b_n\\end{pmatrix}=
        a_1\\cdot b_1+a_2\\cdot b_2+\\dots+a_n\\cdot b_n.

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`

    Returns
    -------
    dot : float
        Dot product value
    """
    # map is noticeably faster than a generator expression for the three
    # element vectors this is almost always called with
    return sum(map(mul, vec_a, vec_b))


def length(vec):
    """Calculate the length of a vector
    :math:`\\boldsymbol{a}\\in\\mathbb{R}^n`

    .. math::

        \\text{length}(\\boldsymbol{a})=|\\boldsymbol{a}|
        =\\sqrt{\\boldsymbol{a}\cdot\\boldsymbol{a}}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    length : float
        Vector length
    """
    return math.sqrt(dot_product(vec, vec))


def vector(pos_a, pos_b):
    """Calculate the vector between to two positions
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{vec}(\\boldsymbol{a},\\boldsymbol{b})
        =\\begin{pmatrix}b_1-a_1\\\\\\vdots\\\\b_n-a_n\\end{pmatrix}

    Parameters
    ----------
    pos_a : list
        First position :math:`\\boldsymbol{a}`
    pos_b : list
        Second position :math:`\\boldsymbol{b}`

    Returns
    -------
    vector : list
        Bond vector
    """
    # Check dimensions
    if not len(pos_a) == len(pos_b):
        print("Vector: Wrong dimensions...")
        return

    # Calculate vector
    return [pos_b[i]-pos_a[i] for i in range(len(pos_a))]


def unit(vec):
    """Transform a vector :math:`\\boldsymbol{a}\\in\\mathbb{R}^n` into a
    unit vector

    .. math::

        \\text{unit}(\\boldsymbol{a})
        =\\frac{\\boldsymbol{a}}{|\\boldsymbol{a}|}

    Parameters
    ----------
    vec : list
        Vector a

    Returns
    -------
    vec : list
        Vector
    """
    vec_length = length(vec)

    return [x/vec_length if not vec_length == 0 else x for x in vec]


def cross_product(vec_a, vec_b):
    """Calculate the cross product of two three-dimensional vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^3`

    .. math::

        \\text{cross}(\\boldsymbol{a},\\boldsymbol{b})=\\begin{pmatrix}
        a_2\\cdot b_3-a_3\\cdot b_2\\\\
        a_3\\cdot b_1-a_1\\cdot b_4\\\\
        a_1\\cdot b_2-a_2\\cdot b_1
        \\end{pmatrix}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`

    Returns
    -------
    vec : list
        Cross product vector
    """
    vec = []
    vec.append(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
    vec.append(vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2])
    vec.append(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

    return vec


def angle(vec_a, vec_b, is_deg=True):
    """Calculate the angle between two vectors
    :math:`\\boldsymbol{a},\\boldsymbol{b}\\in\\mathbb{R}^n`

    .. math::

        \\text{angle}=\\cos^{-1}\\frac{\\boldsymbol{a}\cdot\\boldsymbol{b}}
        {|\\boldsymbol{a}||\\boldsymbol{a}|}

    Parameters
    ----------
    vec_a : list
        First vector :math:`\\boldsymbol{a}`
    vec_b : list
        Second vector :math:`\\boldsymbol{b}`
    is_deg : bool, optional
        True if the output should be in degree

    Returns
    -------
    angle : float
        Angle
    """
    # This is the hot path of the angle sampling, so the three dot products
    # are taken here rather than through dot_product and length, which costs
    # three function calls per angle. The sum(map(mul)) form is kept because
    # it is what dot_product does, down to the summation order
    dot = sum(map(mul, vec_a, vec_b))
    sq_a = sum(map(mul, vec_a, vec_a))
    sq_b = sum(map(mul, vec_b, vec_b))

    # Rounding can push the cosine of parallel vectors marginally outside
    # [-1, 1], which acos rejects, so clamp it to the valid domain
    cos = dot/(math.sqrt(sq_a)*math.sqrt(sq_b))
    angle = math.acos(max(-1.0, min(1.0, cos)))

    return angle*180/math.pi if is_deg else angle
