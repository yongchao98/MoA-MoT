import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * (n choose m) * Gamma(3/2+n)/Gamma(3/2+n-m)
    """
    if n < 0:
        return 0
    
    total_sum = 0
    # Pre-calculate Gamma(3/2 + n) as it's a common term
    try:
        log_gamma_n_plus_3_2 = math.lgamma(n + 1.5)
    except ValueError:
        print(f"Error calculating lgamma({n + 1.5})")
        return float('nan')

    for m in range(n + 1):
        # Calculate combinations (n choose m)
        try:
            comb = math.comb(n, m)
        except ValueError:
            # Should not happen for m <= n
            continue
            
        # Calculate the Gamma function ratio
        # log(Gamma(z1)/Gamma(z2)) = lgamma(z1) - lgamma(z2)
        try:
            log_gamma_ratio = log_gamma_n_plus_3_2 - math.lgamma(n - m + 1.5)
            gamma_ratio = math.exp(log_gamma_ratio)
        except ValueError:
            # This may happen for large arguments if not handled by lgamma
            continue

        # Add the term to the sum
        term = ((-1)**m) * comb * gamma_ratio
        total_sum += term
        
    return total_sum

def main():
    print("Comparing the sum |S_n| with f(n) = n!")
    print("-" * 40)
    print(f"{'n':<5}{'S_n':<20}{'n!':<20}{'|S_n| / n!':<20}")
    print("-" * 40)

    for n in range(11):
        sn = calculate_sum(n)
        nf = math.factorial(n)
        
        # The ratio |S_n|/n! should be bounded by a constant C
        # This is |L_n^(1/2)(1)|
        ratio = abs(sn) / nf if nf > 0 else abs(sn)

        print(f"{n:<5}{sn:<20.6f}{float(nf):<20.1e}{ratio:<20.6f}")

if __name__ == '__main__':
    main()