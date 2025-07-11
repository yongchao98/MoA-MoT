from fractions import Fraction

def solve_optimal_simpson_combination():
    """
    Calculates the constants C, n, m for the error term of the optimal
    linear combination of Simpson's 1/3 and 3/8 rules.
    """
    
    # h is the width of one of the 6 subintervals of [a, b]
    # So, h = (b-a)/6

    # Error series for a single application of the base rules:
    # E(h_base) = C4 * h_base^5 * f^(4) + C6 * h_base^7 * f^(6) + ...

    # For Simpson's 1/3 rule on an interval of width 2h (so h_base = h):
    # Error E_1/3 = -h^5/90 * f^(4) - h^7/1890 * f^(6)
    c4_13 = Fraction(-1, 90)
    c6_13 = Fraction(-1, 1890)

    # For Simpson's 3/8 rule on an interval of width 3h (so h_base = h):
    # Error E_3/8 = -3h^5/80 * f^(4) - 45h^7/448 * f^(6)
    c4_38 = Fraction(-3, 80)
    c6_38 = Fraction(-45, 448)
    
    # The composite 1/3 rule is applied 3 times.
    # The composite 3/8 rule is applied 2 times.
    # Total error for composite rules in terms of h = (b-a)/6:
    # E_comp = N_apps * (C4 * h^5 * f^(4) + C6 * h^7 * f^(6))
    
    # Composite 1/3 rule (3 applications)
    comp_c4_13 = 3 * c4_13
    comp_c6_13 = 3 * c6_13
    
    # Composite 3/8 rule (2 applications)
    comp_c4_38 = 2 * c4_38
    comp_c6_38 = 2 * c6_38

    # We want to find alpha and beta such that:
    # 1. alpha + beta = 1
    # 2. The f^(4) error term is cancelled: alpha*comp_c4_13 + beta*comp_c4_38 = 0
    #
    # From (2), alpha*comp_c4_13 + (1-alpha)*comp_c4_38 = 0
    # alpha * (comp_c4_13 - comp_c4_38) = -comp_c4_38
    # alpha = -comp_c4_38 / (comp_c4_13 - comp_c4_38)
    
    alpha = -comp_c4_38 / (comp_c4_13 - comp_c4_38)
    beta = 1 - alpha

    # The new error term is the sum of the f^(6) terms
    # E_new = (alpha * comp_c6_13 + beta * comp_c6_38) * h^7 * f^(6)
    
    error_coeff_h7 = alpha * comp_c6_13 + beta * comp_c6_38
    
    # The error is E = error_coeff_h7 * h^7 * f^(6)(xi)
    # Substitute h = (b-a)/6
    # E = error_coeff_h7 * ((b-a)/6)^7 * f^(6)(xi)
    # E = (error_coeff_h7 / 6^7) * (b-a)^7 * f^(6)(xi)
    
    n = 7
    m = 6
    C = error_coeff_h7 / (6**n)

    print("The optimal combination is S_new = alpha * S_1/3 + beta * S_3/8")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}")
    print("\nThe error term for this new rule is of the form C * (b-a)^n * f^(m)(xi).")
    print(f"The constants are:")
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"C = {C.numerator}/{C.denominator}")
    
    print("\nThe final error equation is:")
    print(f"Error = ({C.numerator}/{C.denominator}) * (b-a)^{n} * f^({m})(xi)")

    # Return the final tuple for the answer block
    return C, n, m

if __name__ == '__main__':
    C_val, n_val, m_val = solve_optimal_simpson_combination()
    # The final answer format requires the fraction C to be represented as a string
    # and the whole tuple enclosed in <<< >>>
    final_answer = f"<<<({C_val.numerator}/{C_val.denominator}, {n_val}, {m_val})>>>"
    # print(final_answer) # This would be captured for the final answer block

solve_optimal_simpson_combination()