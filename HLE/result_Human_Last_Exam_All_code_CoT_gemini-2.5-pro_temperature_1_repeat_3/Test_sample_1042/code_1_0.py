def demonstrate_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    using a numerical example.

    We aim to identify E(Y^a | A, L) for a=0 and L=1.
    This means we need to find:
    1. E(Y^{a=0} | A=0, L=1)
    2. E(Y^{a=0} | A=1, L=1)
    """

    # --- Assumed Known Quantities ---
    # These values would be computed from data or are given as identifiable.
    # Let's assume a specific value for the confounder, L=1.

    # 1. Given as identifiable (e.g., via some instrumental variable method)
    E_Ya_given_L = 15.0  # E(Y^{a=0} | L=1)

    # 2. Computable from the distribution of observed data (A, Y, L)
    # Probabilities of treatment A conditional on L=1
    P_A1_given_L1 = 0.6
    P_A0_given_L1 = 1 - P_A1_given_L1

    # Observed outcomes conditional on treatment A and L=1
    E_Y_given_A0_L1 = 10.0  # E(Y | A=0, L=1)

    print("Demonstrating identification of E(Y^a | A, L) for a=0 and L=1\n")

    # --- Case 1: Identify E(Y^{a=0} | A=0, L=1) ---
    # This is identified by the consistency assumption.
    E_Y0_given_A0_L1 = E_Y_given_A0_L1
    print("Case 1: Counterfactual treatment matches observed treatment (a=0, A=0)")
    print(f"E(Y^(a=0) | A=0, L=1) is identified by E(Y | A=0, L=1)")
    print(f"E(Y^(a=0) | A=0, L=1) = {E_Y0_given_A0_L1}\n")


    # --- Case 2: Identify E(Y^{a=0} | A=1, L=1) ---
    # This is identified by algebraic manipulation using the law of total expectation.
    # E(Y^a|L) = E(Y^a|A=a,L)P(A=a|L) + E(Y^a|A!=a,L)P(A!=a|L)
    # Solving for E(Y^a|A!=a,L):
    # E(Y^a|A!=a,L) = [ E(Y^a|L) - E(Y^a|A=a,L)P(A=a|L) ] / P(A!=a|L)
    #
    # Plugging in a=0, L=1:
    # E(Y^0|A=1,L=1) = [ E(Y^0|L=1) - E(Y^0|A=0,L=1)P(A=0|L=1) ] / P(A=1|L=1)
    
    numerator = E_Ya_given_L - E_Y0_given_A0_L1 * P_A0_given_L1
    denominator = P_A1_given_L1
    E_Y0_given_A1_L1 = numerator / denominator
    
    print("Case 2: Counterfactual treatment differs from observed treatment (a=0, A=1)")
    print("We use the formula: E(Y^0|A=1,L=1) = [E(Y^0|L=1) - E(Y^0|A=0,L=1) * P(A=0|L=1)] / P(A=1|L=1)")
    print("Plugging in the numbers:")
    print(f"E(Y^0|A=1,L=1) = [ {E_Ya_given_L} - {E_Y0_given_A0_L1} * {P_A0_given_L1} ] / {P_A1_given_L1}")
    print(f"E(Y^0|A=1,L=1) = {numerator} / {denominator}")
    print(f"E(Y^0|A=1,L=1) = {E_Y0_given_A1_L1:.4f}")

if __name__ == '__main__':
    demonstrate_identification()