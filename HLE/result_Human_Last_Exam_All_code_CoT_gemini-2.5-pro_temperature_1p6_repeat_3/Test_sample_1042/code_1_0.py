import pandas as pd

def demonstrate_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L) with a
    numerical example, based on the provided assumptions.

    Assumptions for this numerical example:
    - A is binary {0, 1}, L is binary {0, 1}.
    - We want to identify E(Y^a | A, L) for a=1, i.e., E(Y^1 | A, L).
    
    The script will calculate the four quantities:
    1. E(Y^1 | A=1, L=0)
    2. E(Y^1 | A=1, L=1)
    3. E(Y^1 | A=0, L=0)
    4. E(Y^1 | A=0, L=1)
    """

    print("--- Numerical Demonstration for Identification of E(Y^a | A,L) for a=1 ---\n")

    # --- Step 1: Define quantities assumed to be identifiable ---

    # These would be estimated from the observed data distribution P(Y, A, L)
    # E(Y | A=a, L=l)
    E_Y_given_A1_L0 = 10.0
    E_Y_given_A1_L1 = 15.0
    # P(A=a | L=l)
    P_A1_given_L0 = 0.3
    P_A1_given_L1 = 0.6
    
    # These are derived from the above probabilities
    P_A0_given_L0 = 1 - P_A1_given_L0
    P_A0_given_L1 = 1 - P_A1_given_L1

    # This is the crucial quantity assumed to be identified, perhaps via an
    # instrumental variable or other advanced method.
    # E(Y^a | L=l)
    E_Y1_given_L0 = 9.0
    E_Y1_given_L1 = 16.0

    print("Identifiable Quantities (Given/Assumed):")
    print(f"  E(Y | A=1, L=0) = {E_Y_given_A1_L0}")
    print(f"  E(Y | A=1, L=1) = {E_Y_given_A1_L1}")
    print(f"  P(A=1 | L=0) = {P_A1_given_L0}")
    print(f"  P(A=1 | L=1) = {P_A1_given_L1}")
    print(f"  E(Y^1 | L=0) = {E_Y1_given_L0}")
    print(f"  E(Y^1 | L=1) = {E_Y1_given_L1}")
    print("-" * 20)

    # --- Step 2: Identification for the A=a case ---
    # By consistency, E(Y^1 | A=1, L) = E(Y | A=1, L)
    
    ident_E_Y1_A1_L0 = E_Y_given_A1_L0
    ident_E_Y1_A1_L1 = E_Y_given_A1_L1
    
    print("\nIdentification for A=a (i.e., A=1):")
    print("Based on the consistency axiom: E(Y^1 | A=1, L) = E(Y | A=1, L)")
    print(f"  E(Y^1 | A=1, L=0) = E(Y | A=1, L=0) = {ident_E_Y1_A1_L0}")
    print(f"  E(Y^1 | A=1, L=1) = E(Y | A=1, L=1) = {ident_E_Y1_A1_L1}")
    print("-" * 20)

    # --- Step 3: Identification for the A!=a case ---
    # Using the formula: E(Y^a | A!=a, L) = [E(Y^a|L) - E(Y|A=a,L)P(A=a|L)] / P(A!=a|L)
    
    # For L=0
    numerator_L0 = E_Y1_given_L0 - E_Y_given_A1_L0 * P_A1_given_L0
    ident_E_Y1_A0_L0 = numerator_L0 / P_A0_given_L0
    
    # For L=1
    numerator_L1 = E_Y1_given_L1 - E_Y_given_A1_L1 * P_A1_given_L1
    ident_E_Y1_A0_L1 = numerator_L1 / P_A0_given_L1

    print("\nIdentification for A!=a (i.e., A=0):")
    print("Using formula: E(Y^1 | A=0, L) = [E(Y^1|L) - E(Y|A=1,L)P(A=1|L)] / P(A=0|L)\n")

    print(f"  Calculation for L=0:")
    print(f"    E(Y^1 | A=0, L=0) = (E(Y^1|L=0) - E(Y|A=1,L=0) * P(A=1|L=0)) / P(A=0|L=0)")
    print(f"    E(Y^1 | A=0, L=0) = ({E_Y1_given_L0} - {E_Y_given_A1_L0} * {P_A1_given_L0}) / {P_A0_given_L0:.2f}")
    print(f"    E(Y^1 | A=0, L=0) = ({numerator_L0}) / {P_A0_given_L0:.2f}")
    print(f"    E(Y^1 | A=0, L=0) = {ident_E_Y1_A0_L0:.4f}\n")
    
    print(f"  Calculation for L=1:")
    print(f"    E(Y^1 | A=0, L=1) = (E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)")
    print(f"    E(Y^1 | A=0, L=1) = ({E_Y1_given_L1} - {E_Y_given_A1_L1} * {P_A1_given_L1}) / {P_A0_given_L1:.2f}")
    print(f"    E(Y^1 | A=0, L=1) = ({numerator_L1}) / {P_A0_given_L1:.2f}")
    print(f"    E(Y^1 | A=0, L=1) = {ident_E_Y1_A0_L1:.4f}\n")
    print("-" * 20)

    # --- Step 4: Final results ---
    print("\nSummary of Identified Values for E(Y^1 | A, L):")
    results = {
        ('A=1, L=0'): ident_E_Y1_A1_L0,
        ('A=1, L=1'): ident_E_Y1_A1_L1,
        ('A=0, L=0'): ident_E_Y1_A0_L0,
        ('A=0, L=1'): ident_E_Y1_A0_L1
    }
    for key, val in results.items():
        print(f"  E(Y^1 | {key}) = {val:.4f}")

if __name__ == '__main__':
    demonstrate_identification()