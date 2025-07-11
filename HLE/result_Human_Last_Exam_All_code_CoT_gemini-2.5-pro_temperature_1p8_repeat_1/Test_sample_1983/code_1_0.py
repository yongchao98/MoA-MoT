import math

def solve_difference_equation_problem():
    """
    Solves the problem based on the provided parameters and the formulas
    from the reference paper [1].
    """
    # Parameters given in the problem statement
    # k1 = 10**3000, but we will use its exponent for logarithmic calculations
    # k2 = 10**500, similar to k1
    k1_exp = 3000
    k2_exp = 500
    lambda1 = 0.5
    lambda2 = 0.5
    h_norm = 1000

    # From Theorem 2 in [1], we estimate the limits of the norm of the solution.
    # L_sup = lim_sup ||x_n|| (n -> +inf) = (k2 * |||h|||) / (1 - lambda2)
    # L_inf = lim_inf ||x_n|| (n -> -inf) = (k1 * |||h|||) / (1 - lambda1)

    # We perform calculations using logarithms to handle the large exponents.
    # L_sup = (10^500 * 10^3) / (1 - 0.5) = 10^503 / 0.5 = 2 * 10^503
    # L_inf = (10^3000 * 10^3) / (1 - 0.5) = 10^3003 / 0.5 = 2 * 10^3003

    # log10(L_sup) = log10(2 * 10^503) = log10(2) + 503
    log10_L_sup = math.log10(2) + k2_exp + math.log10(h_norm) - math.log10(1 - lambda2)

    # log10(L_inf) = log10(2 * 10^3003) = log10(2) + 3003
    log10_L_inf = math.log10(2) + k1_exp + math.log10(h_norm) - math.log10(1 - lambda1)

    # The expression to calculate is:
    # E = 100 * log10( (1/3) * L_sup ) + 10 * log10( (1/3) * L_inf )
    # This simplifies to:
    # E = 100 * (log10(L_sup) - log10(3)) + 10 * (log10(L_inf) - log10(3))

    term1 = 100 * (log10_L_sup - math.log10(3))
    term2 = 10 * (log10_L_inf - math.log10(3))
    result = term1 + term2
    
    # We display the final equation with the numbers plugged in.
    # Using log(a/b) = log(a) - log(b), the terms inside the parentheses are:
    # log10(L_sup/3) and log10(L_inf/3)
    # To make it more explicit with numbers:
    # term1 = 100 * (log10(2 * 10^503 / 3)) = 100 * (log10(2/3) + 503)
    # term2 = 10 * (log10(2 * 10^3003 / 3)) = 10 * (log10(2/3) + 3003)

    log10_2_div_3 = math.log10(2.0/3.0)
    L_sup_exp = k2_exp + 3
    L_inf_exp = k1_exp + 3
    
    print("The final expression is:")
    # Using f-string formatting to display the components of the final equation
    final_eq_str = (
        f"100 * (log10(2/3) + {L_sup_exp}) + 10 * (log10(2/3) + {L_inf_exp})"
    )
    print(final_eq_str)

    # Printing the evaluation step by step
    term1_val = 100 * (log10_2_div_3 + L_sup_exp)
    term2_val = 10 * (log10_2_div_3 + L_inf_exp)
    total_val = term1_val + term2_val

    print(f"= ({term1_val}) + ({term2_val})")
    print(f"= {total_val}")


solve_difference_equation_problem()