import numpy as np

def solve():
    """
    This function performs the derivation and prints the final answer.
    The derivation is primarily mathematical and explained in the text above.
    This code serves to formalize the result.
    """
    
    # Let L(n) be the loss term, E[n - tau].
    # We have E[tau] = n - L(n).
    # We need to find theta such that L(n) <= c * n^theta.
    # The derivation shows L(n) is bounded by C * n^(1/2).
    # E[tau] >= n - C * n^(1/2).
    # So, theta = 1/2.
    
    theta_numerator = 1
    theta_denominator = 2
    
    # The question asks for a multiple of 1/8.
    # 1/2 = 4/8.
    
    theta = "1/2"
    
    print("The derivation proceeds as follows:")
    print("1. Express the expectation as E[tau] = n - sum_{j=1}^{n-1} P(S_j >= T).")
    print("2. Bound each probability P(S_j >= T) using Cantelli's inequality (a one-sided Chebyshev inequality).")
    print("3. This requires calculating the mean and variance of S_j.")
    print(f"   E[X_i] = 1/(2n)")
    print(f"   Var(X_i) = 1/(3n^(3/2)) - 1/(4n^2)")
    print(f"   E[S_j] = j/(2n)")
    print(f"   Var(S_j) = j * (1/(3n^(3/2)) - 1/(4n^2))")
    print("4. The bound on P(S_j >= T) is approximately proportional to Var(S_j), which is O(j * n^(-3/2)).")
    print("5. The sum of these probabilities from j=1 to n-1 is approximated by an integral.")
    print("6. The integral is evaluated to be of order O(n^(1/2)).")
    print("7. Therefore, the loss term is bounded by c*n^(1/2).")
    print("8. This means E[tau] >= n - c*n^(1/2), so the largest possible value for theta is 1/2.")
    print("As a multiple of 1/8, theta = 4/8 = 1/2.")
    print("\nFinal Answer:")
    # Final equation format: E[tau] >= n - c*n^(theta)
    # We print the components of this equation.
    print(f"n_term = n")
    print(f"c_term = c")
    print(f"theta = {theta_numerator}/{theta_denominator}")

solve()