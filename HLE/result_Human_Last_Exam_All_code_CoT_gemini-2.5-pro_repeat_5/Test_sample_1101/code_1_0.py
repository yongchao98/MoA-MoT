import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(n+1.5)/Gamma(n-m+1.5)
    """
    total = 0
    # Use log-gamma for numerical stability
    log_gamma_n_plus_1_5 = math.lgamma(n + 1.5)
    for m in range(n + 1):
        # Binomial coefficient C(n, m)
        log_binom_coeff = math.lgamma(n + 1) - math.lgamma(m + 1) - math.lgamma(n - m + 1)
        
        # Gamma ratio term
        log_gamma_n_minus_m_plus_1_5 = math.lgamma(n - m + 1.5)
        log_gamma_ratio = log_gamma_n_plus_1_5 - log_gamma_n_minus_m_plus_1_5
        
        # Combine terms in log scale
        log_term = log_binom_coeff + log_gamma_ratio
        term = math.exp(log_term)
        
        if m % 2 == 1:
            total -= term
        else:
            total += term
            
    return total

def calculate_f(n):
    """
    Calculates the bounding function f(n) = Gamma(n+1.5)
    """
    return math.gamma(n + 1.5)

def main():
    try:
        n_values = [0, 1, 2, 3, 4, 5, 10, 15]
        print("Comparing the sum S_n with the bounding function f(n) = Gamma(n+1.5)")
        print("-" * 60)
        print(f"{'n':>2} | {'S_n':>25} | {'f(n)':>25}")
        print("-" * 60)
        
        for n in n_values:
            s_n = calculate_sum(n)
            f_n = calculate_f(n)
            print(f"{n:>2} | {s_n:>25.4e} | {f_n:>25.4e}")
            
        print("-" * 60)
        
        n_final = 10
        s_n_final = calculate_sum(n_final)
        f_n_final = calculate_f(n_final)
        
        # The inequality is S_n <= C * f(n).
        # We demonstrate that the ratio |S_n/f(n)| is bounded.
        # Let's find a possible C for n=15
        s_15 = calculate_sum(15)
        f_15 = calculate_f(15)
        C = abs(s_15/f_15)

        print(f"The analysis suggests f(n) = Gamma(n+3/2).")
        print("The relationship is sum <= C * f(n).")
        print(f"For n = {n_final}, the sum is {s_n_final:.4e}.")
        print(f"The bounding function f({n_final}) is {f_n_final:.4e}.")
        print(f"We need to show there exists a C > 0 such that |S_n| <= C * f(n).")
        print(f"For n=15, the ratio |S_n/f(n)| is approx {C:.4f}, suggesting such a C exists.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
