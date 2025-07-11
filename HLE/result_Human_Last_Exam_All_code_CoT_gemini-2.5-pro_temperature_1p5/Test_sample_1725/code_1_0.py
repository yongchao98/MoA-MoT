import numpy as np

def solve_crystal_problem():
    """
    Solves for the two requested parameters of the crystal model.
    A. The dimension of the fibers of the bundle pi.
    B. The number of coefficients specifying the energy functional E.
    """

    # Part A: What is the dimension of pi's fibers?
    # The problem describes C as a connection on the tangent bundle T(R^3) that
    # encodes local crystal structure. In the context of continuum mechanics of
    # crystals, this is most standardly interpreted as the distortion tensor field, C_ij.
    # This field measures the local deformation and rotation of the crystal lattice.
    # At each point in space, C_ij is a 3x3 matrix. The space of 3x3 matrices
    # has a dimension of 3 * 3 = 9. This is the dimension of the fiber of the bundle pi,
    # representing the local degrees of freedom of the crystal structure.
    dimension_of_fiber = 9

    # Part B: How many coefficients specify E?
    # The energy density E must satisfy several conditions:
    # 1. Induces linear dynamics -> E must be quadratic in the relevant fields.
    # 2. Detects dislocations (related to curl of C) -> E must depend on derivatives of C.
    # 3. Homogeneous polynomial of least degree -> E is quadratic in the first derivative of C.
    # 4. Invariant under cubic crystal symmetry -> Invariant under the O_h point group.
    #
    # Therefore, E is a quadratic form in the rank-3 tensor T_kij = partial_k(C_ij).
    # E = L_ijklmn * T_kij * T_lmn
    # The number of coefficients is the number of independent components of the
    # tensor L, which must be invariant under the O_h group transformations.
    #
    # This number can be found using group representation theory. It is equal to the
    # number of scalar (A_1g) invariants that can be formed from the tensor product
    # of the representation of T with itself. This number is given by summing the
    # squares of the multiplicities of the irreducible representations (irreps) in the
    # decomposition of the representation of T_kij.
    #
    # The representation of T_kij transforms as a rank-3 tensor. The index 'k' from the
    # derivative transforms as a vector (the T_1u irrep of O_h). The indices 'ij' from
    # the distortion tensor also transform as a tensor product of vectors (T_1u (x) T_1u).
    # Thus, the full representation for T_kij is W = T_1u (x) T_1u (x) T_1u.
    # Since T_1u is an odd ('u') representation, the triple product W is also odd,
    # meaning it can only be composed of 'u' irreps of O_h.

    # We use the character table of the rotational subgroup O (order 24) to find
    # the multiplicities using the character projection formula.
    # The classes of O are: E, 8*C3, 3*C2(C4^2), 6*C4, 6*C2'
    class_sizes = np.array([1, 8, 3, 6, 6])
    group_order = np.sum(class_sizes)  # 24

    # Characters of the 'u' (odd parity) irreps of O_h for the classes in O.
    char_table_u_reps = {
        'A1u': np.array([1,  1,  1, -1, -1]),
        'A2u': np.array([1,  1,  1,  1,  1]),
        'Eu':  np.array([2, -1,  2,  0,  0]),
        'T1u': np.array([3,  0, -1,  1, -1]),
        'T2u': np.array([3,  0, -1, -1,  1]),
    }

    # The character of our representation W is the cube of the vector rep T1u.
    chi_T1u = char_table_u_reps['T1u']
    chi_W = chi_T1u**3

    # The number of quadratic invariants is sum(m_R^2), where m_R is the
    # multiplicity of irrep R in the decomposition of W.
    # m_R = (1/|G|) * sum_over_classes( class_size * chi_R^* * chi_W )
    multiplicities_sq = []
    for chi_R in char_table_u_reps.values():
        # Characters are real, so chi_R^* (conjugate) is just chi_R.
        m_R_numerator = np.sum(class_sizes * chi_R * chi_W)
        m_R = int(round(m_R_numerator / group_order))
        multiplicities_sq.append(m_R**2)

    num_energy_coeffs = sum(multiplicities_sq)

    # Print the final result in the format "A B"
    print(f"{dimension_of_fiber} {num_energy_coeffs}")

solve_crystal_problem()