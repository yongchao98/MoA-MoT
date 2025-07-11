import math

def check_unboundedness_statement(s, p):
    """
    Checks the statement from question (a) for given s and p.

    The statement is: If p > 2(1 + 3s) / (1 + s), then J_t becomes
    unbounded from below as t -> +inf.
    """

    print(f"--- Checking for s = {s}, p = {p} ---")

    # 1. Check if the condition in the question is met
    condition_p_value = 2 * (1 + 3 * s) / (1 + s)
    is_condition_met = p > condition_p_value

    print(f"Condition from question: p > 2*(1+3*s)/(1+s)")
    print(f"Substituting values: {p} > 2*(1 + 3*{s}) / (1 + {s})")
    print(f"Result: {p} > {condition_p_value:.4f}")
    print(f"Is the condition met? {is_condition_met}\n")

    if not is_condition_met:
        print("The premise of the statement is not met for these values.")
        return

    # 2. Determine the actual behavior of the functional J_t
    # The term with the largest exponent of t determines the behavior for t -> inf.
    
    # Exponent of the positive kinetic term
    exp_kin = max(2 * s, 2)
    print("For J_t to go to +inf, the kinetic term must dominate.")
    print(f"The exponent of the dominant kinetic term is max(2*s, 2) = max({2*s}, 2) = {exp_kin}")

    # Exponent of the negative p-term
    exp_pot = (p - 2) * (1 + s) / 2
    print("For J_t to go to -inf, a potential term must dominate.")
    print(f"The exponent of the p-term is (p-2)*(1+s)/2 = ({p}-2)*(1+{s})/2 = {exp_pot}\n")

    # 3. Conclude based on the exponents
    print("Comparing the exponents:")
    if exp_pot > exp_kin:
        print(f"Result: {exp_pot:.4f} > {exp_kin:.4f}. The potential term dominates.")
        print("Conclusion: J_t becomes unbounded from below (goes to -inf). The statement holds for this case.")
    else:
        print(f"Result: {exp_pot:.4f} <= {exp_kin:.4f}. The kinetic term dominates or is equal.")
        print("Conclusion: J_t is bounded below or goes to +inf. The statement FAILS for this case.")


# Demonstrate the counterexample found in the analysis
# We choose s < 1, e.g., s = 0.5.
# For s=0.5, the condition becomes p > 10/3 ~ 3.33. Let's pick p = 4.
s_counterexample = 0.5
p_counterexample = 4

check_unboundedness_statement(s_counterexample, p_counterexample)
