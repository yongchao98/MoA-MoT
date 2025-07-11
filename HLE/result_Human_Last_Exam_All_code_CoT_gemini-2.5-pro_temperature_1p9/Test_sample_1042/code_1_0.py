import numpy as np
import pandas as pd

def demonstrate_identification():
    """
    Demonstrates the identification of E(Y^a | A, L).
    """
    # 1. Define and simulate a Structural Causal Model (SCM)
    # U is an unmeasured binary confounder
    # L is a measured binary confounder
    # A is a binary treatment, caused by L and U
    # Y is the outcome, caused by A, L, and U
    np.random.seed(42)
    n_samples = 1_000_000
    U = np.random.binomial(1, 0.5, n_samples)
    L = np.random.binomial(1, 0.4, n_samples)
    # A is more likely if U=1 or L=1
    prob_A = 1 / (1 + np.exp(-(2 * U + 1.5 * L - 2)))
    A = np.random.binomial(1, prob_A, n_samples)
    # Y is a linear function of A, L, U + noise
    Y = 2 * A + 3 * L + 4 * U + np.random.normal(0, 1, n_samples)

    # The potential outcomes Y^a are what Y would be if A were set to a.
    Y0 = 2 * 0 + 3 * L + 4 * U + np.random.normal(0, 1, n_samples) # Y if A=0
    Y1 = 2 * 1 + 3 * L + 4 * U + np.random.normal(0, 1, n_samples) # Y if A=1

    df = pd.DataFrame({'U': U, 'L': L, 'A': A, 'Y': Y, 'Y0': Y0, 'Y1': Y1})

    print("--- Identification Demonstration ---\n")
    # Loop over interventions a=0 and a=1
    for a in [0, 1]:
        print(f"--- Intervention: A = {a} ---")
        # 2. Calculate the TRUE value of E(Y^a | A, L) using the full SCM data (incl. U)
        # This is what we want to identify, but cannot compute directly without U.
        true_values = df.groupby(['A', 'L'])[[f'Y{a}']].mean().rename(columns={f'Y{a}': 'True E[Y^a|A,L]'})

        # 3. Emulate the researcher: use only observed data (A, L, Y) and the identifiable quantity E(Y^a|L)
        
        # As per the problem, E(Y^a | L) is assumed to be identifiable. We calculate it
        # from the full simulation to serve as the "oracle-provided" quantity.
        oracle_E_Ya_L = df.groupby('L')[[f'Y{a}']].mean()

        # The researcher can calculate these from observed data (df[['A', 'L', 'Y']])
        p_A_given_L = df.groupby('L')['A'].value_counts(normalize=True).rename('p(A|L)').reset_index()
        E_Y_given_A_L = df.groupby(['A', 'L'])['Y'].mean().rename('E(Y|A,L)').reset_index()

        # 4. Apply the identification formulas
        results = []
        for l_val in [0, 1]:
            for a_prime in [0, 1]: # a_prime is the conditioning value of A
                # Get necessary components from observed data calculations
                E_Ya_L_val = oracle_E_Ya_L.loc[l_val, f'Y{a}']
                p_a_L = p_A_given_L.query(f"L == {l_val} and A == {a}").iloc[0]['p(A|L)']
                p_not_a_L = 1 - p_a_L
                
                print(f"\nCalculating E[Y^{a} | A={a_prime}, L={l_val}]:")
                
                # Case 1: a_prime == a (Identification by consistency rule)
                if a_prime == a:
                    E_Y_a_L = E_Y_given_A_L.query(f"A == {a} and L == {l_val}").iloc[0]['E(Y|A,L)']
                    ident_val = E_Y_a_L
                    print(f"  Using consistency E(Y^a | A=a, L) = E(Y | A=a, L)")
                    print(f"  E[Y | A={a}, L={l_val}] = {ident_val:.4f}")

                # Case 2: a_prime != a (Identification using the derived formula)
                else: 
                    E_Y_a_L = E_Y_given_A_L.query(f"A == {a} and L == {l_val}").iloc[0]['E(Y|A,L)']
                    ident_val = (E_Ya_L_val - p_a_L * E_Y_a_L) / p_not_a_L
                    print(f"  Using formula: (E[Y^a|L] - P(A=a|L) * E[Y|A=a,L]) / P(A!=a|L)")
                    print(f"  E[Y^{a} | A={a_prime}, L={l_val}] = ({E_Ya_L_val:.4f} - {p_a_L:.4f} * {E_Y_a_L:.4f}) / {p_not_a_L:.4f}")
                    print(f"  Result = {ident_val:.4f}")

                true_val = true_values.loc[(a_prime, l_val), 'True E[Y^a|A,L]']
                results.append({
                    'Intervention (a)': a,
                    'Condition A': a_prime,
                    'Condition L': l_val,
                    'Identified Value': ident_val,
                    'True Value': true_val
                })
        
        results_df = pd.DataFrame(results)
        print("\n--- Comparison of Identified vs. True values ---")
        print(results_df.to_string(index=False))
        print("-" * 50)


if __name__ == '__main__':
    demonstrate_identification()
