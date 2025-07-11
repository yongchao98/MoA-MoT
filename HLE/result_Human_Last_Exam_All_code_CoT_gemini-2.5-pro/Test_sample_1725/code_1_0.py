def solve_crystal_problem():
    """
    This function calculates the dimension of the fiber and the number of energy coefficients
    based on the physics and mathematics of the provided continuum model.
    """

    # Part A: Dimension of pi's fibers
    # The field C is a connection on the tangent bundle of R^3.
    # Its components at a point are given by a tensor C^i_{jk}.
    # The indices i, j, k each run from 1 to 3.
    # The dimension is the total number of components.
    num_indices = 3
    rank = 3
    fiber_dimension = num_indices ** rank

    # Part B: Number of coefficients specifying E
    # E is a quadratic form in dC, E ~ (dC)^2.
    # The number of coefficients is the number of independent quadratic invariants
    # of the tensor dC under the cubic symmetry group O_h.
    # This requires a group theory calculation based on the decomposition of the
    # representation of dC into irreducible representations of the cubic group.
    #
    # The representation of dC (a tensor D^i_{j,kl}) under SO(3) is 1 tensor 1 tensor 1.
    # This decomposes as V_SO(3) = 0 + 3(1) + 2(2) + 3.
    # We then decompose this into irreps of the cubic group O:
    # V_O = A_1 + A_2 + 2E + 4T_1 + 3T_2
    # The number of invariants is sum over R of m_R * (m_R + 1) / 2,
    # where m_R is the multiplicity of the irrep R.
    
    # Multiplicities of irreps A1, A2, E, T1, T2
    multiplicities = {
        'A1': 1,
        'A2': 1,
        'E': 2,
        'T1': 4,
        'T2': 3,
    }

    num_coeffs = 0
    for m in multiplicities.values():
        num_coeffs += m * (m + 1) // 2
    
    energy_coefficients = num_coeffs

    print(f"{fiber_dimension} {energy_coefficients}")

solve_crystal_problem()