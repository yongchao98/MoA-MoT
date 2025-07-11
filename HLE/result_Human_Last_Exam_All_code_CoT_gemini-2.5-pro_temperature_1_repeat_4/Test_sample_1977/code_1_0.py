import math

def calculate_T_norm_for_even_n(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for a given even integer n.
    """
    if n % 2 != 0 or n < 0:
        print(f"Error: n must be a non-negative even integer. Received n = {n}.")
        return

    # The derived formula for the 1-norm of T is:
    # 3 * 2**(n+1) + 4**(n+1) - (2**(n+3) * (2**n + 1)) / (1 + 3**n)

    # Calculate each part of the formula
    term1_base = 2
    term1_exp = n + 1
    term1_val = 3 * (term1_base**term1_exp)

    term2_base = 4
    term2_exp = n + 1
    term2_val = term2_base**term2_exp

    num_term1_base = 2
    num_term1_exp = n + 3
    num_term1_val = num_term1_base**num_term1_exp

    num_term2_base = 2
    num_term2_exp = n
    num_term2_val = num_term2_base**num_term2_exp + 1
    
    numerator_val = num_term1_val * num_term2_val

    den_term_base = 3
    den_term_exp = n
    denominator_val = 1 + den_term_base**den_term_exp
    
    term3_val = numerator_val / denominator_val

    # Final result
    result = term1_val + term2_val - term3_val

    # Print the step-by-step calculation
    print(f"For n = {n}:")
    print("The 1-norm is calculated using the formula:")
    print("3 * 2**(n+1) + 4**(n+1) - (2**(n+3) * (2**n + 1)) / (1 + 3**n)\n")
    
    print("Plugging in n = " + str(n) + ":")
    # Step 1: Show the formula with n substituted
    print(f"3 * {term1_base}**({n}+{term1_exp-n}) + {term2_base}**({n}+{term2_exp-n}) - ({num_term1_base}**({n}+{num_term1_exp-n}) * ({num_term2_base}**{n} + 1)) / (1 + {den_term_base}**{n})")
    
    # Step 2: Show the exponents evaluated
    print(f"= 3 * {term1_base}**{term1_exp} + {term2_base}**{term2_exp} - ({num_term1_base}**{num_term1_exp} * ({num_term2_base**num_term2_exp} + 1)) / (1 + {den_term_base**den_term_exp})")

    # Step 3: Show the powers evaluated
    print(f"= 3 * {term1_base**term1_exp} + {term2_base**term2_exp} - ({num_term1_val} * {num_term2_val}) / {denominator_val}")
    
    # Step 4: Show the multiplications
    print(f"= {term1_val} + {term2_val} - {numerator_val} / {denominator_val}")

    # Step 5: Show the division
    print(f"= {term1_val} + {term2_val} - {term3_val}")
    
    # Step 6: Show the addition
    print(f"= {term1_val + term2_val} - {term3_val}")

    # Step 7: Final result
    print(f"= {result}")

# Example calculation for n=2
calculate_T_norm_for_even_n(2)