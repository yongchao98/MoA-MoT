import numpy as np

def solve():
    """
    Solves the two parts of the problem.
    Part A is based on physical interpretation.
    Part B is calculated using group theory.
    """

    # Part A: Dimension of the fiber
    # The connection C is an R^3-valued 1-form on R^3.
    # It has components C^a_i where a, i in {1, 2, 3}.
    # This gives 3x3 = 9 components.
    dim_pi_fibers = 9

    # Part B: Number of coefficients for E
    # This is determined by the number of quadratic invariants of the field T ~ dC
    # under the cubic symmetry group O_h. The field T transforms as T_1u x T_1u.
    # We use character theory to find the number of irreps in its decomposition.

    # Character table for O_h group.
    # Classes: E, 8C3, 6C2, 6C4, 3C2', i, 8S6, 6sd, 6S4, 3sh
    # The order of classes can vary, but the inner product is invariant.
    # We use a standard ordering and corresponding characters.
    # Class orders (number of elements in each class)
    class_orders = np.array([1, 8, 6, 6, 3, 1, 8, 6, 6, 3])
    group_order = np.sum(class_orders)

    # Irreducible representations (irreps) of O_h
    irreps = {
        'A1g': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        'A2g': np.array([1, 1, -1, -1, 1, 1, 1, -1, -1, 1]),
        'Eg':  np.array([2, -1, 0, 0, 2, 2, -1, 0, 0, 2]),
        'T1g': np.array([3, 0, -1, 1, -1, 3, 0, -1, 1, -1]),
        'T2g': np.array([3, 0, 1, -1, -1, 3, 0, 1, -1, -1]),
        'A1u': np.array([1, 1, 1, 1, 1, -1, -1, -1, -1, -1]),
        'A2u': np.array([1, 1, -1, -1, 1, -1, -1, 1, 1, -1]),
        'Eu':  np.array([2, -1, 0, 0, 2, -2, 1, 0, 0, -2]),
        'T1u': np.array([3, 0, -1, 1, -1, -3, 0, 1, -1, 1]),
        'T2u': np.array([3, 0, 1, -1, -1, -3, 0, -1, 1, 1]),
    }

    # The field T transforms as the tensor product T1u x T1u.
    # The character of a tensor product is the product of the characters.
    char_T1u = irreps['T1u']
    char_product = char_T1u * char_T1u

    # Decompose the product representation into irreps by computing the inner product
    # of its character with each irrep's character.
    # multiplicity(i) = (1/|G|) * sum_g(chi_prod(g) * chi_i(g)^*)
    multiplicities = {}
    for name, char_irrep in irreps.items():
        # Characters are real for O_h
        inner_product = np.sum(class_orders * char_product * char_irrep)
        multiplicity = int(round(inner_product / group_order))
        if multiplicity > 0:
            multiplicities[name] = multiplicity

    # The number of coefficients is the number of irreps in the decomposition,
    # assuming multiplicity one for each. If multiplicities were > 1, the
    # number of coefficients would be sum(m_i^2).
    # Let's check the decomposition:
    # T1u x T1u = A1g + Eg + T1g + T2g
    # All multiplicities are 1.
    num_coeffs = 0
    for m in multiplicities.values():
        num_coeffs += m**2
        
    # The final answer is a pair of numbers.
    print(f"{dim_pi_fibers} {num_coeffs}")

solve()