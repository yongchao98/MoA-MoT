import numpy as np

def solve():
    """
    Solves the problem by finding the eigenvalues of a specific matrix
    and calculating the infimum and supremum based on them.
    """
    # Step 1: Find the real root of the polynomial x^3 + x^2 + x - 1 = 0
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    xi = None
    for r in roots:
        if np.isreal(r):
            xi = r.real
            break

    if xi is None:
        print("Error: Real root not found.")
        return

    # Step 2: Construct the matrix S based on expressions in xi
    # Using polynomial relations to simplify the entries:
    # 1/xi = xi^2+xi+1
    # 1/xi^2 = 4*xi^2+4*xi
    s11 = 2
    s12 = -1 - xi
    s13 = -xi**2 - 1
    s21 = s12
    s22 = 2 * (xi**2 + xi + 1)
    s23 = -(xi**2 + xi + 2)
    s31 = s13
    s32 = s23
    s33 = 8 * (xi**2 + xi)
    
    S = np.array([
        [s11, s12, s13],
        [s21, s22, s23],
        [s31, s32, s33]
    ])

    # Step 3: Compute the eigenvalues of S
    eigenvalues = np.linalg.eigvalsh(S)
    
    # Eigenvalues are sorted in increasing order by np.linalg.eigvalsh
    # S is a rank-2 matrix, so one eigenvalue should be close to zero.
    lambda_min_nonzero = eigenvalues[1]
    lambda_max = eigenvalues[2]

    # Step 4: Calculate infimum and supremum
    inf_val = 2 / lambda_max
    sup_val = 2 / lambda_min_nonzero

    print(f"The real root is xi = {xi}")
    print("\nThe matrix S is:")
    print(S)
    print(f"\nThe eigenvalues of S are: {eigenvalues}")
    print(f"\ninf_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = 2 / {lambda_max:.8f} = {inf_val:.8f}")
    print(f"sup_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = 2 / {lambda_min_nonzero:.8f} = {sup_val:.8f}")

solve()
<<<
inf = 0.23751219
sup = 0.49138195
>>>