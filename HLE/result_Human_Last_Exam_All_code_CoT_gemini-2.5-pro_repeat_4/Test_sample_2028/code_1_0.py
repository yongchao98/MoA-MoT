import numpy as np
from itertools import combinations

def solve_vest_subset_sum():
    """
    Solves an example instance of the VEST problem for commuting matrices.
    This demonstrates the calculation described in the problem, for which an FPT
    algorithm exists as explained in part (a).
    """
    # (a) Example: Commuting matrices (diagonal), S=I
    # Let n=2, m=4, k=2.
    n = 2
    m = 4
    k = 2

    # v is a rational vector
    v = np.array([1, 2])

    # S is the identity matrix
    S = np.identity(n)

    # T_i are commuting rational matrices. We use diagonal matrices for simplicity.
    T = [
        np.diag([1, 2]),
        np.diag([3, 1]),
        np.diag([2, 2]),
        np.diag([1, 3])
    ]

    print(f"Solving VEST for k={k} with m={m} commuting matrices.")

    indices = list(range(m))
    subsets_indices = combinations(indices, k)
    
    total_sum = 0
    terms = []

    for I_indices in subsets_indices:
        # Product of matrices in the subset I
        P_I = np.identity(n)
        for i in I_indices:
            P_I = P_I @ T[i]

        # Since matrices are diagonal, P_I^T = P_I
        P_I_T = P_I.T

        # Calculate the term v^T * P_I * S * P_I^T * v
        term = v.T @ P_I @ S @ P_I_T @ v
        terms.append(term)
        total_sum += term
        
    # Output the equation with each number
    equation_str = " + ".join(map(str, terms))
    print(f"The sum is: {equation_str} = {total_sum}")

if __name__ == '__main__':
    solve_vest_subset_sum()