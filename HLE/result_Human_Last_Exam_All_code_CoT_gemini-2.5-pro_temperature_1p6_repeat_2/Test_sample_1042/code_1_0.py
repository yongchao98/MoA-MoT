import pandas as pd

def solve_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L) using a
    numerical example.

    It calculates the values of E(Y^a | A, L) for a=1 and for a binary A and L,
    based on quantities assumed to be identifiable from data or by assumption.
    """

    # --- Assumed Identifiable Quantities ---
    # These values would be estimated from observational data or given.
    # We want to find E(Y^1 | A, L)
    
    # Let's assume A is binary (0, 1) and L is binary (0, 1).
    
    # 1. From observational data (Y, A, L)
    # Conditional expectations of observed outcome Y (identifiable via regression)
    E_Y_given_A1_L0 = 10.0
    E_Y_given_A1_L1 = 15.0
    
    # Propensity scores (identifiable from data)
    P_A1_given_L0 = 0.4
    P_A0_given_L0 = 1 - P_A1_given_L0
    
    P_A1_given_L1 = 0.7
    P_A0_given_L1 = 1 - P_A1_given_L1
    
    # 2. By the problem's primary assumption
    # Average counterfactual outcome conditional on L
    E_Y1_given_L0 = 12.0
    E_Y1_given_L1 = 16.0

    # --- Identification Calculation ---
    
    # Case 1: a=A (e.g., Y^1 for those with A=1)
    # E(Y^1 | A=1, L=l) = E(Y | A=1, L=l) by consistency
    
    E_Y1_given_A1_L0 = E_Y_given_A1_L0
    print("Identification using consistency rule:")
    print(f"E(Y^1 | A=1, L=0) = E(Y | A=1, L=0) = {E_Y1_given_A1_L0}")
    
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    print(f"E(Y^1 | A=1, L=1) = E(Y | A=1, L=1) = {E_Y1_given_A1_L1}")
    print("-" * 30)

    # Case 2: a!=A (e.g., Y^1 for those with A=0)
    # E(Y^1 | A=0,L) = [E(Y^1|L) - E(Y^1|A=1,L)P(A=1|L)] / P(A=0|L)

    # For L=0
    numerator_L0 = E_Y1_given_L0 - E_Y1_given_A1_L0 * P_A1_given_L0
    E_Y1_given_A0_L0 = numerator_L0 / P_A0_given_L0
    print("Identification using algebraic rearrangement:")
    print(f"E(Y^1 | A=0, L=0) = (E(Y^1|L=0) - E(Y^1|A=1,L=0) * P(A=1|L=0)) / P(A=0|L=0)")
    print(f"                 = ({E_Y1_given_L0} - {E_Y1_given_A1_L0} * {P_A1_given_L0}) / {P_A0_given_L0}")
    print(f"                 = ({E_Y1_given_L0 - E_Y1_given_A1_L0 * P_A1_given_L0}) / {P_A0_given_L0}")
    print(f"                 = {E_Y1_given_A0_L0:.4f}")

    # For L=1
    numerator_L1 = E_Y1_given_L1 - E_Y1_given_A1_L1 * P_A1_given_L1
    E_Y1_given_A0_L1 = numerator_L1 / P_A0_given_L1
    print(f"\nE(Y^1 | A=0, L=1) = (E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)) / P(A=0|L=1)")
    print(f"                 = ({E_Y1_given_L1} - {E_Y1_given_A1_L1} * {P_A1_given_L1}) / {P_A0_given_L1}")
    print(f"                 = ({E_Y1_given_L1 - E_Y1_given_A1_L1 * P_A1_given_L1}) / {P_A0_given_L1}")
    print(f"                 = {E_Y1_given_A0_L1:.4f}")

solve_identification()