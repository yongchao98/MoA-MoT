import numpy as np

def calculate_ngsc(eigenvalues, condition_name):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC) from a set of eigenvalues.
    NGSC is the normalized Shannon entropy of the normalized eigenvalues.
    This function also prints the components of the calculation.
    """
    print(f"--- Calculating NGSC for '{condition_name}' condition ---")
    
    # Ensure eigenvalues are positive to avoid errors with log(0).
    eigenvalues = np.array(eigenvalues)
    eigenvalues[eigenvalues <= 0] = 1e-12 
    
    # 1. Normalize eigenvalues to unit sum (p_i = lambda_i / sum(lambda))
    total_variance = np.sum(eigenvalues)
    normalized_eigenvalues = eigenvalues / total_variance
    num_components = len(eigenvalues)
    
    print(f"Number of components (N): {num_components}")
    
    # 2. Calculate Shannon entropy (H = -sum(p_i * log(p_i)))
    # We use the natural logarithm (base e).
    entropy = -np.sum(normalized_eigenvalues * np.log(normalized_eigenvalues))
    
    # 3. Calculate maximum possible entropy (H_max = log(N))
    max_entropy = np.log(num_components)
    
    # 4. Normalize the entropy to get NGSC (NGSC = H / H_max)
    ngsc_value = entropy / max_entropy
    
    print("\nEquation for NGSC: H / H_max")
    print(f"Entropy (H) = {entropy:.4f}")
    print(f"Max Entropy (H_max = log({num_components})) = {max_entropy:.4f}")
    print(f"Final Equation: NGSC = {entropy:.4f} / {max_entropy:.4f}")
    print(f"Resulting NGSC for '{condition_name}': {ngsc_value:.4f}\n")
    return ngsc_value

# --- Simulation based on Figure 3a (right) ---
# We simulate eigenvalues for 50 components, matching the figure's x-axis.

N_COMPONENTS = 50

# Case 1: Pre-Psilocybin (e.g., MTP/No Drug). Eigenvalues decay quickly.
# This corresponds to the steep blue/grey curves, indicating high synchrony.
components_pre = np.arange(1, N_COMPONENTS + 1)
eigenvalues_pre = 1 / (components_pre**2)

# Case 2: Post-Psilocybin. Eigenvalues decay slowly.
# This corresponds to the shallower red curve, indicating higher spatial complexity.
components_post = np.arange(1, N_COMPONENTS + 1)
eigenvalues_post = 1 / components_post

# Calculate and display NGSC for both simulated conditions.
ngsc_pre = calculate_ngsc(eigenvalues_pre, "Pre-Psilocybin (High Synchrony)")
ngsc_post = calculate_ngsc(eigenvalues_post, "Post-Psilocybin (High Desynchrony)")

# --- Conclusion from Simulation ---
print("--- Analysis ---")
print("The simulation shows that a more even distribution of variance across components (slower eigenvalue decay, as in the 'Post-Psilocybin' case) yields a higher NGSC value.")
print(f"Simulated Pre-Psilocybin NGSC: {ngsc_pre:.4f}")
print(f"Simulated Post-Psilocybin NGSC: {ngsc_post:.4f}")
print("This demonstrates a significant increase in NGSC, which is consistent with the findings in Figure 3b and supports answer choice G.")
