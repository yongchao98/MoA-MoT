import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n) / Gamma(1.5+n-m)
    """
    if n < 0:
        return 0
    
    # The term Gamma(1.5+n) can be pulled out of the sum.
    try:
        gamma_n_plus_1_5 = math.gamma(1.5 + n)
    except ValueError:
        # Handle cases where the gamma function is not defined, though not expected for n in N
        return float('nan')

    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n, m)
        comb = math.comb(n, m)
        
        # Calculate the Gamma term in the denominator
        try:
            gamma_denom = math.gamma(1.5 + n - m)
        except ValueError:
            # This case can happen if 1.5 + n - m is zero or a negative integer
            # which is not possible for m <= n
            gamma_denom = float('inf')

        term = (-1)**m * comb / gamma_denom
        total_sum += term
        
    return total_sum * gamma_n_plus_1_5

def analyze_complexity(max_n=15):
    """
    Analyzes the complexity of the sum by comparing it to Gamma(n+1.5).
    """
    print(f"{'n':>2} | {'S_n':>20} | {'f(n)=Gamma(n+1.5)':>20} | {'|S_n| / f(n)':>15}")
    print("-" * 65)
    
    for n in range(max_n + 1):
        sn = calculate_sum(n)
        
        # Candidate function f(n) = Gamma(n+1.5)
        fn = math.gamma(n + 1.5)
        
        ratio = abs(sn) / fn if fn != 0 else float('inf')
        
        print(f"{n:>2} | {sn:>20.6f} | {fn:>20.6f} | {ratio:>15.6f}")

# Run the analysis
analyze_complexity()
