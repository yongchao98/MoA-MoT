import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n) / Gamma(3/2+n-m)
    for a given integer n.
    """
    if n < 0:
        print("n must be a non-negative integer.")
        return

    total_sum = 0
    
    # Pre-calculate Gamma(3/2 + n)
    try:
        gamma_n_plus_3_2 = math.gamma(1.5 + n)
    except ValueError:
        print(f"Error calculating gamma function for n={n}. n might be too large.")
        return

    for m in range(n + 1):
        # Calculate binomial coefficient C(n, m)
        comb = math.comb(n, m)
        
        # Calculate falling factorial as Gamma(3/2+n)/Gamma(3/2+n-m)
        try:
            # For m=0, the denominator gamma is the same as numerator
            if m == 0:
                term_val = 1.0
            else:
                gamma_n_plus_3_2_minus_m = math.gamma(1.5 + n - m)
                term_val = gamma_n_plus_3_2 / gamma_n_plus_3_2_minus_m
        except ValueError:
            # This can happen if 1.5 + n - m is a non-positive integer.
            # However, for m <= n, 1.5 + n - m is always positive.
            # This handles potential overflow for large n.
            print(f"Error calculating gamma function for n={n}, m={m}.")
            continue

        # Add the term to the sum
        term = ((-1)**m) * comb * term_val
        total_sum += term
        
    print(f"The sum S_n for n = {n} is:")
    print(total_sum)
    # The sum is bounded by C*f(n) where f(n) = n!
    # So we print the final equation showing the sum and the bounding function f(n)=n!
    print("\nThe final equation is:")
    # Re-calculate parts to print the equation term-by-term for the first few terms
    # to avoid excessive output length.
    max_terms_to_print = min(n + 1, 5) 
    equation_str = []
    for m in range(max_terms_to_print):
        sign = "-" if m % 2 != 0 else ""
        if m > 0:
            sign = " - " if m % 2 != 0 else " + "
        else:
            sign = ""

        comb = math.comb(n, m)
        gamma_ratio = math.gamma(1.5 + n) / math.gamma(1.5 + n - m)
        
        equation_str.append(f"{sign}{comb}*{gamma_ratio:.4f}")
    
    if n + 1 > max_terms_to_print:
        equation_str.append("...")
        
    final_sum_str = f" = {total_sum}"
    
    print(" ".join(equation_str) + final_sum_str)


# Example usage for n=4
n = 4
calculate_sum(n)