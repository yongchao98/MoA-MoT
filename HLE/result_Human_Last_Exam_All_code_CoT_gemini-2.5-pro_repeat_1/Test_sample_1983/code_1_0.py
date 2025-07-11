import math

def solve_dichotomy_problem():
    """
    Solves the problem based on the theory of discrete dichotomies.
    """
    
    # Step 1: Define constants from the problem
    k1_exp = 3000
    k2_exp = 500
    lambda1 = 0.5
    lambda2 = 0.5
    h_norm = 1000

    # Step 2: Calculate the effective limit values L_+ and L'_-
    # L_+ = (k1 / (1 - lambda1)) * h_norm
    # L_plus_val_mantissa = (1 / (1 - lambda1)) * h_norm = (1 / 0.5) * 1000 = 2000 = 2 * 10^3
    # L_plus_exp = k1_exp + 3 = 3003
    # L_+ = 2 * 10^3003
    
    # L'_- = (k2 * lambda2 / (1 - lambda2)) * h_norm
    # L_minus_val_mantissa = (lambda2 / (1 - lambda2)) * h_norm = (0.5 / 0.5) * 1000 = 1000 = 1 * 10^3
    # L_minus_exp = k2_exp + 3 = 503
    # L'_- = 1 * 10^503

    log10_L_plus_mantissa = math.log10(2)
    log10_L_plus_exp = 3003
    
    log10_L_minus_mantissa = math.log10(1) # This is 0
    log10_L_minus_exp = 503

    log10_L_plus = log10_L_plus_mantissa + log10_L_plus_exp
    log10_L_minus = log10_L_minus_mantissa + log10_L_minus_exp
    
    # Step 3: Calculate the final expression
    # E = 100 * log10(L_+ / 3) + 10 * log10(L'_- / 3)
    # E = 100 * (log10(L_+) - log10(3)) + 10 * (log10(L'_-) - log10(3))
    
    log10_3 = math.log10(3)
    
    term1 = 100 * (log10_L_plus - log10_3)
    term2 = 10 * (log10_L_minus - log10_3)
    
    result = term1 + term2

    # Step 4: Print the final equation and the result
    final_equation = f"100*log10( (1/3) * (2 * 10^{log10_L_plus_exp}) ) + 10*log10( (1/3) * (1 * 10^{log10_L_minus_exp}) )"
    print("The final equation with the computed values is:")
    print(final_equation)
    print("\nThe numerical result is:")
    print(result)

solve_dichotomy_problem()