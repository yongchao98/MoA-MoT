import sys

def solve():
    """
    Solves the two parts of the crystal continuum model problem.

    Part A: What is the dimension of pi's fibers?

    The physical state of the crystal at a point is described by the connection C, which in
    continuum mechanics corresponds to the distortion tensor, C_ij. This is a general
    3x3 real matrix, relating the deformed state to a reference state. The space of
    such matrices at each point forms the fiber of the bundle pi. The dimension of the
    space of 3x3 real matrices is 3 * 3 = 9.
    """
    dimension_of_fiber = 9

    """
    Part B: How many coefficients specify E?

    The energy density E is a homogeneous polynomial of the least degree that is
    invariant under rigid transformations and the crystal's cubic symmetry. The
    degrees are defined as deg(C) = 1 and deg(partial) = 1. The term with the
    least degree is quadratic in C, with degree 2 (e.g., C_ij * C_kl). This is
    the elastic energy. A term involving the dislocation density dC (proportional to partial C)
    would have a higher degree (e.g., (dC)^2 has degree 4).

    So, E must be a quadratic function of the components of C. We need to find the number of
    independent coefficients in this quadratic form that are invariant under cubic symmetry.

    We decompose C into its symmetric (strain) part, epsilon, and its anti-symmetric
    (rotation) part, omega. The energy is a sum of contributions from each:
    E(C) = E_symm(epsilon) + E_anti(omega).

    1. E_symm(epsilon): This is the standard elastic energy. For a cubic crystal, there are 3
       independent elastic constants (e.g., c11, c12, c44).
    2. E_anti(omega): The anti-symmetric tensor omega has 3 independent components and transforms
       as an axial vector. The only quadratic form invariant under cubic symmetry is the
       squared magnitude, which requires 1 coefficient.

    Total coefficients = 3 (from symmetric part) + 1 (from anti-symmetric part) = 4.
    """
    number_of_coefficients = 4
    
    # The final answer requires printing each number separated by a space.
    print(f"{dimension_of_fiber} {number_of_coefficients}")

solve()
<<<9 4>>>