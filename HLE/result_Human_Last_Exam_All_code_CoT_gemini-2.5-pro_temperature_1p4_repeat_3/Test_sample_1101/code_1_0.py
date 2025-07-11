import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n) / Gamma(1.5+n-m)
    """
    if n < 0:
        return 0
    
    # Use log-gamma for numerical stability with larger n
    log_gamma_n_plus_1_5 = math.lgamma(1.5 + n)
    
    total_sum = 0
    for m in range(n + 1):
        # The term can be written as [n+1/2]_m which is a product.
        # This is more stable than computing gamma functions of large numbers.
        term = 1.0
        for i in range(m):
            term *= (n + 0.5 - i)
            
        # Binomial coefficient
        binom = math.comb(n, m)
        
        # Add or subtract from total
        if m % 2 == 1:
            total_sum -= binom * term
        else:
            total_sum += binom * term
            
    return total_sum

def f(n):
    """
    Calculates the bounding function f(n) = n! / n^(1/4)
    """
    if n == 0:
        # Define f(0) as 1 to avoid division by zero and match S_0
        return 1.0
    # Use log-factorial (lgamma(n+1)) for stability
    log_f_n = math.lgamma(n + 1) - 0.25 * math.log(n)
    return math.exp(log_f_n)

def main():
    """
    Main function to perform the calculation and print the results.
    """
    print("This script finds the function f(n) with the lowest complexity such that")
    print("the sum |S_n| <= C*f(n) for a constant C > 0.")
    print("\nBased on asymptotic analysis, we propose f(n) = n! / n^(1/4).")
    print("We will now verify this numerically by computing the ratio |S_n/f(n)|,\nwhich should be bounded by C.\n")

    running_max_C = 0.0
    
    print(f"{'n':<5}{'S_n':<20}{'f(n)':<25}{'|S_n / f(n)| (<= C)':<25}")
    print("-" * 75)
    
    for n in range(15):
        try:
            sum_val = calculate_sum(n)
            f_val = f(n)
            
            if abs(f_val) < 1e-9: # Avoid division by zero
                ratio = float('inf')
            else:
                ratio = abs(sum_val / f_val)
            
            if ratio != float('inf') and ratio > running_max_C:
                running_max_C = ratio

            # The final equation is |S_n| <= C * f(n).
            # We print the components of this inequality.
            print(f"{n:<5}{sum_val:<20.4f}{f_val:<25.4e}{ratio:<25.4f}")

        except (ValueError, OverflowError) as e:
            print(f"{n:<5}Error calculating: {e}")
            break
            
    print("-" * 75)
    print(f"\nThe numerical results suggest the ratio is bounded.")
    print(f"A possible value for C is the maximum ratio found: {running_max_C:.4f}")

if __name__ == '__main__':
    main()
