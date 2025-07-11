import numpy as np

def calculate_ngsc(eigenvalues):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC).

    Args:
        eigenvalues (list or np.ndarray): A list or array of non-negative eigenvalues from PCA.

    Returns:
        float: The NGSC value.
    """
    eigenvalues = np.array(eigenvalues, dtype=float)
    
    # The total variance is the sum of eigenvalues.
    total_variance = np.sum(eigenvalues)
    if total_variance == 0:
        return 0.0

    # The normalized eigenvalues are the proportion of variance for each component.
    # p_i = eigenvalue_i / total_variance
    normalized_eigenvalues = eigenvalues / total_variance
    
    # We only consider components with non-zero variance to avoid log(0).
    non_zero_normalized_eigenvalues = normalized_eigenvalues[normalized_eigenvalues > 0]
    
    num_components = len(non_zero_normalized_eigenvalues)
    if num_components <= 1:
        # If variance is in 0 or 1 component, entropy is 0.
        return 0.0

    # Shannon entropy S = - sum(p_i * log(p_i))
    # We use natural log (np.log) as is standard for S_max = log(N).
    entropy = -np.sum(non_zero_normalized_eigenvalues * np.log(non_zero_normalized_eigenvalues))

    # Maximum possible entropy is log(N), where N is the number of components.
    max_entropy = np.log(num_components)

    # NGSC is the entropy normalized by the maximum possible entropy.
    ngsc = entropy / max_entropy
    return ngsc

# --- Evaluate Theoretical Claims ---
print("--- Evaluating Theoretical Answer Choices ---")

# For M: NGSC = 0 interpretation.
# This represents perfect synchrony (all variance in one component).
eigs_sync = [100, 0, 0, 0]
ngsc_sync = calculate_ngsc(eigs_sync)
print(f"Option M Check: For eigenvalues = {eigs_sync} (perfect synchrony), the calculation is NGSC = 0 / log(1) which is 0.")
print(f"Calculated NGSC = {ngsc_sync:.4f}. This is a state of minimum complexity (perfect synchrony), not signal independence. Option M is incorrect.")
print("-" * 30)

# For D: Even distribution of variance.
# This represents maximum complexity/desynchronization.
eigs_even = [25, 25, 25, 25]
norm_eigs_even = [eig / sum(eigs_even) for eig in eigs_even]
entropy_even = -sum([p * np.log(p) for p in norm_eigs_even])
max_entropy_even = np.log(len(eigs_even))
ngsc_even = calculate_ngsc(eigs_even)
print(f"Option D Check: For eigenvalues = {eigs_even} (even distribution), the entropy is -4 * (0.25 * log(0.25)) = {entropy_even:.4f}.")
print(f"The max entropy is log(4) = {max_entropy_even:.4f}. Thus NGSC = {entropy_even:.4f} / {max_entropy_even:.4f} = {ngsc_even:.4f}.")
print("An NGSC of 1 represents maximum complexity, which implies low, not high, functional connectivity. Option D is incorrect.")
print("-" * 30)

# For E, H, N: Range of NGSC.
print(f"Option E, H, N Check: The NGSC is a normalized entropy value.")
print(f"Its minimum value (perfect synchrony) is {ngsc_sync:.1f}. Its maximum value (perfect desynchronization) is {ngsc_even:.1f}.")
print("Therefore, NGSC is always in the range [0, 1]. It cannot be 2 (E), negative (H), or -1 (N). These options are incorrect.")
print("-" * 30)

# --- Evaluate Figure-Based Claims ---
print("--- Evaluating Figure-Based Answer Choices ---")
print("1. After invalidating the theoretical options (D,E,H,M,N) and other incorrect claims (B,F,I,J,L) through direct observation, we are left with G and K.")
print("2. Option G: 'Post-psilocybin whole-brain NGSC shows a significant increase...'. This is strongly supported by Fig 3b, where red dots (psilocybin) are consistently and clearly higher than grey/blue dots for every participant. This represents the main finding of the plot.")
print("3. Option K: 'Participant 4 has more evenly distributed data variance... under each psilocybin... scan than any other... scan'. This is also correct. In Fig 3b, all of P4's red dots (NGSC > ~0.78) are visually higher than all of P4's grey and blue dots (NGSC < ~0.73).")
print("\nConclusion: Since both statements G and K are correct and directly supported by the provided figure, the correct answer choice is C.")