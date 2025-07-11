import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n)/Gamma(3/2+n-m)
    and prints the full equation and the result.
    
    The term Gamma(a)/Gamma(a-m) is calculated as the product (a-1)*(a-2)*...*(a-m).
    """
    if not isinstance(n, int) or n < 0:
        print("Please provide a non-negative integer for n.")
        return

    # Using math.lgamma for precision with large numbers
    log_gamma_n_plus_3_2 = math.lgamma(n + 1.5)

    total_sum = 0.0
    equation_parts = []

    for m in range(n + 1):
        # Calculate C(n, m)
        log_comb = math.lgamma(n + 1) - math.lgamma(m + 1) - math.lgamma(n - m + 1)
        comb = round(math.exp(log_comb))

        # Calculate Gamma ratio
        log_gamma_ratio = log_gamma_n_plus_3_2 - math.lgamma(n + 1.5 - m)
        gamma_ratio = math.exp(log_gamma_ratio)

        term = ((-1)**m) * comb * gamma_ratio
        total_sum += term
        
        # Build the equation string
        # We explicitly show the numbers in the equation
        gamma_ratio_val = 1.0
        for i in range(m):
            gamma_ratio_val *= (n + 0.5 - i)
            
        sign = "-" if m % 2 != 0 else "+"
        if m == 0:
            equation_parts.append(f"{comb} * {gamma_ratio_val:.4f}")
        else:
            equation_parts.append(f" {sign} {comb} * {gamma_ratio_val:.4f}")
    
    print(f"For n = {n}, the equation is:")
    # Join parts; handle the initial plus sign for display
    full_equation = "".join(equation_parts).lstrip(" +")
    print(f"S_{n} = {full_equation}")
    print(f"Result: S_{n} = {total_sum}")

# Example execution for n=4
calculate_sum(4)

# Example execution for a slightly larger n
calculate_sum(5)