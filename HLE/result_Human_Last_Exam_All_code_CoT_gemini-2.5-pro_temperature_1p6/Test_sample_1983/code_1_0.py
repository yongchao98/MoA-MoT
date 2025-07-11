import math

def solve_dichotomy_problem():
    """
    Solves the problem based on the provided parameters for a difference equation
    admitting a discrete dichotomy.
    """
    
    # Step 1: Interpret the constants
    k1_exp = 3000
    k2_exp = 500
    h_norm = 1000
    
    # Interpreting lambda_2 = 0.5 * lambda_1 = 0.5 as:
    # 0.5 * lambda_1 = 0.5  => lambda_1 = 1.0
    # lambda_2 = 0.5
    lambda_1 = 1.0
    lambda_2 = 0.5

    # Step 2: Determine Asymptotic Norms
    # The asymptotic norm for n -> +infinity
    # lim_sup_norm = (k1 / lambda_1) * h_norm
    # lim_sup_norm = (10^3000 / 1.0) * 1000 = 10^3003
    log10_lim_sup_norm = k1_exp + math.log10(h_norm / lambda_1)

    # The asymptotic norm for n -> -infinity
    # lim_inf_norm = (k2 / lambda_2) * h_norm
    # lim_inf_norm = (10^500 / 0.5) * 1000 = 2 * 10^503
    log10_lim_inf_norm = k2_exp + math.log10(h_norm / lambda_2)
    
    # Step 3 & 4: Evaluate the Final Expression
    # The expression is: 100 * log10(lim_sup_norm / 3) + 10 * log10(lim_inf_norm / 3)
    #
    # Term 1: 100 * (log10(lim_sup_norm) - log10(3))
    # Term 2: 10 * (log10(lim_inf_norm) - log10(3))
    
    log10_3 = math.log10(3)
    
    term1_val = 100 * (log10_lim_sup_norm - log10_3)
    term2_val = 10 * (log10_lim_inf_norm - log10_3)
    
    result = term1_val + term2_val

    # The problem asks to output the final equation with numbers
    # Final form: A + B * log10(2) - C * log10(3)
    # 100 * (3003 - log10(3)) + 10 * (log10(2*10^503) - log10(3))
    # = 300300 - 100*log10(3) + 10 * (log10(2) + 503 - log10(3))
    # = 300300 - 100*log10(3) + 10*log10(2) + 5030 - 10*log10(3)
    # = 305330 + 10*log10(2) - 110*log10(3)
    
    val_A = 305330
    val_B = 10
    val_C = 110
    
    print("The final calculation is based on the equation:")
    print(f"{val_A} + {val_B} * log10(2) - {val_C} * log10(3)")
    
    final_result_from_formula = val_A + val_B * math.log10(2) - val_C * math.log10(3)
    
    print("\nResult:")
    print(final_result_from_formula)

solve_dichotomy_problem()