def solve_problem():
    """
    This function calculates the maximum number of roots for the given problem
    and prints the step-by-step derivation and the final answer.
    """

    # Part (a): Find the expression for the maximum number of roots in terms of t.

    # The expression for the maximum number of roots is t * (t - 1) / 2.
    # Here's the derivation explained in print statements.

    print("Derivation for the maximum number of roots:")
    print("1. The number of roots of R_t in ]0, 1[ is given by the degree of the Wronskian of polynomials g_i.")
    print("   Number of Roots = sum(deg(g_i)) - t*(t-1)/2")
    
    print("\n2. We choose the canonical set of exponents k_i and l_i to be {0, 1, ..., t-1}.")
    t_str = "t"
    sum_of_exponents_formula = f"{t_str}*({t_str}-1)/2"
    
    print(f"   sum(k_i for i=1..t) = sum(j for j=0..t-1) = {sum_of_exponents_formula}")
    print(f"   sum(l_i for i=1..t) = sum(j for j=0..t-1) = {sum_of_exponents_formula}")
    
    sum_of_degrees_formula = f"2 * ({sum_of_exponents_formula}) = {t_str}*({t_str}-1)"
    print(f"   sum(deg(g_i)) = sum(k_i) + sum(l_i) = {sum_of_degrees_formula}")

    max_roots_formula = f"{sum_of_degrees_formula} - {sum_of_exponents_formula} = {sum_of_exponents_formula}"
    print(f"\n3. Substituting this into the formula for the number of roots:")
    print(f"   Max number of roots = {max_roots_formula}")
    
    final_expression_a = "t*(t-1)/2"

    # Part (b): Calculate this maximum number for t = 5.
    t = 5
    sum_k_l = t * (t - 1)
    t_term = t * (t - 1) // 2
    max_roots_val = sum_k_l - t_term

    print("\nFor t = 5:")
    print(f"Sum of degrees = {t}*({t}-1) = {sum_k_l}")
    print(f"The second term = {t}*({t}-1)/2 = {t_term}")
    print(f"Maximum number of roots = {sum_k_l} - {t_term} = {max_roots_val}")
    
    final_expression_b = str(max_roots_val)
    
    # Final Answer Formatting
    final_answer = f"(a) {final_expression_a}; (b) {final_expression_b}"
    print("\nFinal Answer:")
    print(f"<<<{final_answer}>>>")

solve_problem()