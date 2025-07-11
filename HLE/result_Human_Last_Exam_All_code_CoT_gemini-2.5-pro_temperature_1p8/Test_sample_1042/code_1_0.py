import pandas as pd

def identify_e_ya_given_a_l():
    """
    Demonstrates the identification of E(Y^a | A, L) given E(Y^a | L).

    This function uses a sample observed dataset (A, L, Y) and an assumed
    identifiable quantity E(Y^a | L) to calculate E(Y^a | A, L) for a=1.
    """
    # 1. Sample Observed Data
    # A: Binary treatment {0, 1}
    # L: Binary measured confounder {0, 1}
    # Y: Continuous outcome
    data = {
        'A': [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1],
        'L': [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        'Y': [5, 6, 7, 6, 15, 17, 8, 9, 20, 22, 21, 23]
    }
    df = pd.DataFrame(data)
    print("--- Observed Data ---")
    print(df)
    print("\n")

    # 2. Premise: We can identify E(Y^a | L)
    # Let's assume a=1. The following values for E(Y^1 | L) are given.
    # (These values would be identified through other means, e.g., an instrumental variable).
    E_Y1_given_L = {
        0: 16.0,  # E(Y^1 | L=0)
        1: 21.0   # E(Y^1 | L=1)
    }
    print("--- Given Identifiable Quantities (for a=1) ---")
    print(f"E(Y^1 | L=0) = {E_Y1_given_L[0]}")
    print(f"E(Y^1 | L=1) = {E_Y1_given_L[1]}")
    print("\n")
    
    print("--- Identification of E(Y^1 | A, L) ---")

    # We want to identify the function E(Y^1 | A, L). We do this for each value of (A,L).
    # We iterate through each stratum of L.
    for l_val in [0, 1]:
        print(f"\n--- Stratum L={l_val} ---")
        
        # Case 1: A = a (i.e., A = 1)
        # We use the consistency rule: E(Y^1 | A=1, L=l) = E(Y | A=1, L=l)
        
        # Calculate E(Y | A=1, L=l) from data
        group_A1_L = df[(df['A'] == 1) & (df['L'] == l_val)]
        E_Y_given_A1_L = group_A1_L['Y'].mean()

        print(f"Case 1 (A=a): Identification by Consistency")
        print(f"E(Y^1 | A=1, L={l_val}) = E(Y | A=1, L={l_val}) = {E_Y_given_A1_L:.2f}")

        # Case 2: A != a (i.e., A = 0)
        # We use the law of total expectation:
        # E(Y^1 | A=0, L=l) = [ E(Y^1 | L=l) - E(Y^1 | A=1, L=l) * P(A=1 | L=l) ] / P(A=0 | L=l)

        # Calculate P(A|L) from data
        stratum_L = df[df['L'] == l_val]
        P_A1_given_L = len(stratum_L[stratum_L['A'] == 1]) / len(stratum_L)
        P_A0_given_L = 1 - P_A1_given_L

        if P_A0_given_L == 0:
            print(f"Cannot identify for A=0 as P(A=0 | L={l_val}) is zero (no data).")
            continue
            
        # Get the necessary values
        # E(Y^1 | L=l) is given
        E_Y1_L_val = E_Y1_given_L[l_val]
        # E(Y^1 | A=1, L=l) was just calculated
        E_Y1_A1_L_val = E_Y_given_A1_L

        # Apply the formula
        numerator = E_Y1_L_val - (E_Y1_A1_L_val * P_A1_given_L)
        E_Y1_A0_L_val = numerator / P_A0_given_L
        
        print(f"\nCase 2 (A!=a): Identification by Law of Total Expectation")
        print("Formula: E(Y^1|A=0,L={l}) = [E(Y^1|L={l}) - E(Y|A=1,L={l})*P(A=1|L={l})] / P(A=0|L={l})".format(l=l_val))
        print(f"Calculation:")
        print(f"E(Y^1 | A=0, L={l_val}) = ({E_Y1_L_val} - {E_Y1_A1_L_val:.2f} * {P_A1_given_L:.2f}) / {P_A0_given_L:.2f}")
        print(f"               = {numerator:.2f} / {P_A0_given_L:.2f} = {E_Y1_A0_L_val:.2f}")


if __name__ == '__main__':
    identify_e_ya_given_a_l()
