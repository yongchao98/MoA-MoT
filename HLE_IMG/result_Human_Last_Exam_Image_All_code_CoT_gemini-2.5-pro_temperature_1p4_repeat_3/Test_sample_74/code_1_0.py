import numpy as np

def calculate_ngsc(eigenvalues, case_name):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC)
    from a list of eigenvalues.

    NGSC is the normalized Shannon entropy of the normalized eigenvalues.
    """
    print(f"--- Calculating NGSC for: {case_name} ---")
    eigs = np.array(eigenvalues, dtype=float)
    
    # 1. Normalize eigenvalues to unit sum
    # This represents the proportion of variance explained by each component.
    # Let's call them p_i.
    norm_eigs = eigs / np.sum(eigs)
    
    # 2. Calculate Shannon Entropy (H)
    # H = -Σ(p_i * log(p_i))
    # We filter out p_i = 0, as 0 * log(0) is defined as 0.
    non_zero_eigs = norm_eigs[norm_eigs > 0]
    entropy_h = -np.sum(non_zero_eigs * np.log(non_zero_eigs))
    
    # 3. Calculate Maximum Possible Entropy (H_max)
    # This occurs when variance is spread completely evenly (all p_i are equal).
    # H_max = log(N), where N is the number of components.
    num_components = len(eigs)
    if num_components <= 1:
        # NGSC is 0 if there's only one component (max synchrony).
        max_entropy_h_max = 0
        ngsc = 0
    else:
        max_entropy_h_max = np.log(num_components)
        ngsc = entropy_h / max_entropy_h_max

    # Print the equation components as requested
    print(f"Eigenvalues representing variance: {eigenvalues}")
    print(f"Number of components (N): {num_components}")
    print(f"Normalized eigenvalues (p_i): {np.round(norm_eigs, 3)}")
    
    # Building the entropy equation string
    entropy_calc_str = " + ".join([f"(-{p:.3f} * log({p:.3f}))" for p in non_zero_eigs])
    print(f"Entropy (H) = {entropy_calc_str} = {entropy_h:.4f}")
    
    # Building the max entropy equation string
    max_entropy_calc_str = f"log(N) = log({num_components})"
    print(f"Max Entropy (H_max) = {max_entropy_calc_str} = {max_entropy_h_max:.4f}")

    # Building the final NGSC equation string
    ngsc_calc_str = f"H / H_max = {entropy_h:.4f} / {max_entropy_h_max:.4f}"
    print(f"Final NGSC = {ngsc_calc_str} = {ngsc:.4f}\n")


# Case 1: High Synchrony (Low Complexity)
# Most variance is concentrated in the first few components.
# This corresponds to a low NGSC.
high_sync_eigs = [90, 5, 2, 1, 1, 1]
calculate_ngsc(high_sync_eigs, "High Synchrony / Low Complexity (e.g., No Drug)")

# Case 2: High Desynchrony (High Complexity)
# Variance is spread more evenly across many components.
# This corresponds to a high NGSC.
high_desync_eigs = [20, 18, 16, 15, 12, 19]
calculate_ngsc(high_desync_eigs, "High Desynchrony / High Complexity (e.g., Psilocybin)")
