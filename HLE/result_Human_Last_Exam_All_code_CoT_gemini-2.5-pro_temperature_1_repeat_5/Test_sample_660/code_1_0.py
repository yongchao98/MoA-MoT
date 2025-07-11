import numpy as np

def solve_auditor_choice():
    """
    This function solves for the set of firm types (theta) that are more truthful
    and choose more lenient auditors.

    The solution is based on a revised model where the lenient auditor has a small
    (epsilon) probability of detecting malpractice, making the problem non-trivial.
    """

    # --- Step 1: Define Model Parameters ---
    # We choose numerical values that satisfy the problem's constraints:
    # - t, p, F are decreasing functions of theta.
    # - The sorting condition p[1] < epsilon * p[0] must hold.

    # t(theta): Tax liability. t(0) > t(1)
    t = {0: 100, 1: 50}
    # p(x): Probability of RS audit. p(0) > p(1)
    p = {0: 0.8, 1: 0.05}
    # F(theta): Penalty. F(0) > F(1)
    F = {0: 1000, 1: 100}
    # Epsilon: Small probability for lenient auditor refusal
    epsilon = 0.1

    # Verify the sorting condition
    if not (p[1] < epsilon * p[0]):
        print("Warning: Parameters do not satisfy the sorting condition p[1] < epsilon * p[0].")
        print(f"p[1] = {p[1]}, epsilon*p[0] = {epsilon * p[0]}")
        return

    # --- Step 2: Define Cost Functions based on the Revised Model ---
    def cost_auditor_0(theta):
        """Expected cost of choosing the lenient auditor (x=0)."""
        return t[theta] + epsilon * (1 - theta) * p[0] * F[theta]

    def cost_auditor_1(theta):
        """Expected cost of choosing the strict auditor (x=1)."""
        return t[theta] + (1 - theta) * p[1] * F[theta]

    # --- Step 3: Determine Optimal Auditor for Each Firm Type ---
    firm_types = [0, 1]
    optimal_auditor = {}
    costs = {}

    print("--- Calculating Expected Costs ---")
    for theta in firm_types:
        costs[(theta, 0)] = cost_auditor_0(theta)
        costs[(theta, 1)] = cost_auditor_1(theta)
        print(f"Firm theta={theta}:")
        print(f"  Cost with Lenient Auditor (x=0): {costs[(theta, 0)]:.2f}")
        print(f"  Cost with Strict Auditor (x=1): {costs[(theta, 1)]:.2f}")

        # Tie-breaking rule: if costs are equal, choose lenient auditor (x=0)
        if costs[(theta, 1)] < costs[(theta, 0)]:
            optimal_auditor[theta] = 1 # Strict
        else:
            optimal_auditor[theta] = 0 # Lenient

    print("\n--- Optimal Auditor Choice ---")
    for theta in firm_types:
        choice = "Strict (x=1)" if optimal_auditor[theta] == 1 else "Lenient (x=0)"
        print(f"Firm theta={theta} chooses: {choice}")


    # --- Step 4: Find the Set of More Truthful Firms Choosing More Lenient Auditors ---
    result_set = set()
    # We compare the more truthful firm (theta=1) to the less truthful one (theta=0)
    theta_truthful = 1
    theta_malpractice = 0

    is_more_truthful = theta_truthful > theta_malpractice
    chooses_more_lenient = optimal_auditor[theta_truthful] < optimal_auditor[theta_malpractice]

    if is_more_truthful and chooses_more_lenient:
        result_set.add(theta_truthful)

    print("\n--- Final Result ---")
    print("The set of values for theta for which companies keeping more truthful accounts choose more lenient auditors is:")
    print(result_set)

solve_auditor_choice()