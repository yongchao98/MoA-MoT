import numpy as np

def identify_all_components(L_val):
    """
    This function demonstrates the identification of E(Y^a | A, L) for a given
    stratum of L, based on quantities available from observational data and
    the key assumption that E(Y^a | L) is identifiable.
    
    The numerical values used are derived from a structural causal model in the
    thought process, but here they are treated as given inputs.
    """

    # --- INPUTS ---
    # The following values are specific to the chosen stratum, L = L_val.
    if L_val == 0:
        # Quantities identifiable from observational data P(Y, A, L)
        p_A1_given_L = 0.5
        E_Y_given_A1_L = 3.193
        E_Y_given_A0_L = 0.807
        
        # Quantities assumed to be identifiable (given by an "oracle")
        E_Y1_given_L = 2.5
        E_Y0_given_L = 1.5
    else:
        # These are sample values for another stratum, say L=1
        p_A1_given_L = 0.842
        E_Y_given_A1_L = 4.698
        E_Y_given_A0_L = 2.447

        E_Y1_given_L = 4.5
        E_Y0_given_L = 3.5

    p_A0_given_L = 1 - p_A1_given_L
    
    print(f"--- Demonstrating Identification for L = {L_val} ---")
    
    # 1. Identification via Consistency Rule
    print("\nPart 1: Identifying E[Y^a | A=a, L] using Consistency")
    E_Y1_A1_L = E_Y_given_A1_L
    E_Y0_A0_L = E_Y_given_A0_L
    print(f"E[Y^a=1 | A=1, L={L_val}] is identified as E[Y | A=1, L={L_val}] = {E_Y1_A1_L}")
    print(f"E[Y^a=0 | A=0, L={L_val}] is identified as E[Y | A=0, L={L_val}] = {E_Y0_A0_L}")

    # 2. Identification via Algebraic Manipulation
    print("\nPart 2: Identifying E[Y^a | A!=a, L] using algebra")

    # Calculation for E[Y^1 | A=0, L]
    numerator_1 = E_Y1_given_L - E_Y_given_A1_L * p_A1_given_L
    E_Y1_A0_L = numerator_1 / p_A0_given_L

    print("\nCalculating E[Y^a=1 | A=0, L=0]:")
    print(f"  Formula: (E[Y^1|L] - E[Y|A=1,L] * P(A=1|L)) / P(A=0|L)")
    print(f"  Substituting values: ({E_Y1_given_L} - {E_Y_given_A1_L} * {p_A1_given_L}) / {p_A0_given_L}")
    print(f"  Numerator becomes: ({E_Y1_given_L} - {E_Y_given_A1_L * p_A1_given_L:.3f}) = {numerator_1:.3f}")
    print(f"  Result: {numerator_1:.3f} / {p_A0_given_L} = {E_Y1_A0_L:.3f}")

    # Calculation for E[Y^0 | A=1, L]
    numerator_0 = E_Y0_given_L - E_Y_given_A0_L * p_A0_given_L
    E_Y0_A1_L = numerator_0 / p_A1_given_L
    
    print("\nCalculating E[Y^a=0 | A=1, L=0]:")
    print(f"  Formula: (E[Y^0|L] - E[Y|A=0,L] * P(A=0|L)) / P(A=1|L)")
    print(f"  Substituting values: ({E_Y0_given_L} - {E_Y_given_A0_L} * {p_A0_given_L}) / {p_A1_given_L}")
    print(f"  Numerator becomes: ({E_Y0_given_L} - {E_Y_given_A0_L * p_A0_given_L:.3f}) = {numerator_0:.3f}")
    print(f"  Result: {numerator_0:.3f} / {p_A1_given_L} = {E_Y0_A1_L:.3f}")

if __name__ == '__main__':
    # Run the demonstration for the L=0 stratum
    identify_all_components(L_val=0)
