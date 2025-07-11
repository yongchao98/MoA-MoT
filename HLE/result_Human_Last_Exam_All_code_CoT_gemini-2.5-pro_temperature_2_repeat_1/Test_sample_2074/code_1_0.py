import numpy as np
from scipy.linalg import eigh

def solve():
    """
    Solves the problem by numerically computing the eigenvalues and the objective function.
    """
    n = 101
    b = 0.5

    # Step 1: Construct the matrix T' = (1-b^2)*C(b).
    # Based on the derivation, C(b) is a symmetric tridiagonal matrix.
    # C_{11} = C_{n,n} = 1/(1-b^2)
    # C_{ii} = (1+b^2)/(1-b^2) for 1 < i < n
    # C_{i,i+1} = -b/(1-b^2)
    # So T' = (1-b^2)C(b) is:
    T_prime = np.zeros((n, n))
    
    # Set diagonal elements
    diag_val = 1 + b**2
    T_prime.fill_diagonal(diag_val)
    T_prime[0, 0] = 1
    T_prime[n - 1, n - 1] = 1

    # Set off-diagonal elements
    off_diag_val = -b
    for i in range(n - 1):
        T_prime[i, i + 1] = off_diag_val
        T_prime[i + 1, i] = off_diag_val

    # Step 2: Compute eigenvalues of C(b).
    # Eigenvalues of T' (gamma) are computed first.
    # eigh returns sorted eigenvalues for symmetric matrices.
    gamma = eigh(T_prime, eigvals_only=True)
    
    # Eigenvalues of C(b) are gamma / (1-b^2)
    lambdas = gamma / (1 - b**2)

    # Step 3: Compute min_k [k*lambda_k + sum_{i=k+1 to n} lambda_i]
    # In python lambdas[k-1] corresponds to math lambda_k.
    
    # Precomputing cumulative sums from the right speeds up the loop.
    sums_from_right = np.cumsum(lambdas[::-1])[::-1]

    s_values = []
    for k in range(1, n + 1):  # k from 1 to 101
        # term1 = k * lambda_k
        term1 = k * lambdas[k - 1]
        
        # term2 = sum of lambda_i for i > k
        term2 = sums_from_right[k] if k < n else 0
        
        s_values.append(term1 + term2)

    min_s = min(s_values)

    # Step 4: Compute ell(1/2)
    ell_half = 101 + min_s
    
    # The eigenvalues of C(b) and C(-b) are identical, so ell(b) = ell(-b)
    ell_neg_half = ell_half
    
    # Step 5: Compute the final result
    final_answer = 6 * (ell_half + ell_neg_half)
    
    print(f"The calculated value for ell(1/2) is: {ell_half}")
    print(f"The calculated value for ell(-1/2) is: {ell_neg_half}")
    print(f"The final result is 6 * ({ell_half} + {ell_neg_half}) = {final_answer}")
    
    # Finally, printing the answer in the requested format
    print(f"\n<<< {final_answer} >>>")

solve()
