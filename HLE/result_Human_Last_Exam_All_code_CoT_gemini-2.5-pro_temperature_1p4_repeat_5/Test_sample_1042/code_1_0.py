import pandas as pd

def solve_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    based on the premises in the problem.

    We assume a binary treatment A in {0, 1} and a binary measured
    confounder L in {0, 1}. We want to identify E(Y^1 | A, L).

    This requires identifying four values:
    1. E(Y^1 | A=0, L=0)
    2. E(Y^1 | A=0, L=1)
    3. E(Y^1 | A=1, L=0)
    4. E(Y^1 | A=1, L=1)
    """

    # --- Step 1: Assume values for identifiable quantities ---
    # These would typically be estimated from a large observational dataset.

    # Conditional probabilities of treatment A given confounder L
    P_A1_given_L0 = 0.4
    P_A0_given_L0 = 1 - P_A1_given_L0

    P_A1_given_L1 = 0.7
    P_A0_given_L1 = 1 - P_A1_given_L1

    # Observable expectations E(Y | A, L)
    # By consistency, E(Y | A=1, L=l) = E(Y^1 | A=1, L=l)
    E_Y_given_A1_L0 = 0.8  # This is E(Y^1 | A=1, L=0)
    E_Y_given_A1_L1 = 0.9  # This is E(Y^1 | A=1, L=1)

    # Per the problem's premise, we assume E(Y^a | L) is identified.
    # Let's assume these values were found using some valid method
    # (e.g., from an instrumental variable analysis).
    E_Y1_given_L0 = 0.70
    E_Y1_given_L1 = 0.85

    print("--- Assumed Identifiable Quantities ---")
    print(f"P(A=1|L=0) = {P_A1_given_L0}")
    print(f"P(A=0|L=0) = {P_A0_given_L0}")
    print(f"P(A=1|L=1) = {P_A1_given_L1}")
    print(f"P(A=0|L=1) = {P_A0_given_L1}")
    print(f"E(Y^1|L=0) = {E_Y1_given_L0} (given as identified)")
    print(f"E(Y^1|L=1) = {E_Y1_given_L1} (given as identified)\n")


    # --- Step 2: Identify E(Y^1 | A, L) ---

    # Case 1: A=1 (Observed treatment matches counterfactual)
    # These are directly identified by consistency.
    E_Y1_given_A1_L0 = E_Y_given_A1_L0
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    
    print("--- Identification Step-by-Step ---")
    print("Finding E(Y^1 | A=1, L=0) and E(Y^1 | A=1, L=1):")
    print("By consistency, E(Y^1 | A=1, L) = E(Y | A=1, L)")
    print(f"Identified E(Y^1 | A=1, L=0) = {E_Y1_given_A1_L0}")
    print(f"Identified E(Y^1 | A=1, L=1) = {E_Y1_given_A1_L1}\n")

    # Case 2: A=0 (Observed treatment differs from counterfactual)
    # We use the formula derived from the law of total expectation.
    # E(Y^1|A=0,L) = [E(Y^1|L) - E(Y^1|A=1,L)P(A=1|L)] / P(A=0|L)
    
    print("Finding E(Y^1 | A=0, L=0):")
    # Calculation for L=0
    numerator_L0 = E_Y1_given_L0 - E_Y1_given_A1_L0 * P_A1_given_L0
    denominator_L0 = P_A0_given_L0
    E_Y1_given_A0_L0 = numerator_L0 / denominator_L0
    
    print("The equation is: E(Y^1|A=0,L=0) = [E(Y^1|L=0) - E(Y^1|A=1,L=0) * P(A=1|L=0)] / P(A=0|L=0)")
    print(f"Plugging in the numbers: ({E_Y1_given_L0} - {E_Y1_given_A1_L0} * {P_A1_given_L0}) / {P_A0_given_L0}")
    print(f"Identified E(Y^1 | A=0, L=0) = {E_Y1_given_A0_L0:.4f}\n")

    # Calculation for L=1
    numerator_L1 = E_Y1_given_L1 - E_Y1_given_A1_L1 * P_A1_given_L1
    denominator_L1 = P_A0_given_L1
    E_Y1_given_A0_L1 = numerator_L1 / denominator_L1
    
    print("Finding E(Y^1 | A=0, L=1):")
    print("The equation is: E(Y^1|A=0,L=1) = [E(Y^1|L=1) - E(Y^1|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)")
    print(f"Plugging in the numbers: ({E_Y1_given_L1} - {E_Y1_given_A1_L1} * {P_A1_given_L1}) / {P_A0_given_L1}")
    print(f"Identified E(Y^1 | A=0, L=1) = {E_Y1_given_A0_L1:.4f}\n")


    # --- Step 3: Display Final Result ---
    # We have now identified all four quantities that constitute E(Y^1|A,L)
    
    results = {
        ('L=0', 'A=0'): E_Y1_given_A0_L0,
        ('L=0', 'A=1'): E_Y1_given_A1_L0,
        ('L=1', 'A=0'): E_Y1_given_A0_L1,
        ('L=1', 'A=1'): E_Y1_given_A1_L1,
    }
    
    results_df = pd.DataFrame(pd.Series(results, name='E(Y^1|A,L)'),
                              index=pd.MultiIndex.from_tuples(results.keys(), names=['L', 'A']))

    print("--- Final Identified Values for E(Y^1|A,L) ---")
    print(results_df)

if __name__ == '__main__':
    solve_identification()