import numpy as np

def calculate_ngsc(eigenvalues):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC).

    Args:
        eigenvalues (list or np.array): A list of eigenvalues from PCA.

    Returns:
        float: The NGSC value.
    """
    # Ensure eigenvalues are positive and convert to numpy array
    eigenvalues = np.array(eigenvalues)
    if np.any(eigenvalues < 0):
        raise ValueError("Eigenvalues cannot be negative.")
    
    # Normalize eigenvalues to unit sum (p_i)
    total_variance = np.sum(eigenvalues)
    if total_variance == 0:
        return 0 # No variance, no complexity
    
    normalized_eigenvalues = eigenvalues / total_variance
    
    # Filter out zero-valued eigenvalues to avoid log(0)
    non_zero_eigenvalues = normalized_eigenvalues[normalized_eigenvalues > 0]
    
    # Number of components
    N = len(eigenvalues)
    if N <= 1:
        return 0 # Entropy is 0 if there's only one component
        
    # Calculate Shannon entropy (H)
    shannon_entropy = -np.sum(non_zero_eigenvalues * np.log(non_zero_eigenvalues))
    
    # Calculate maximum possible entropy (log(N))
    max_entropy = np.log(N)
    
    # Calculate NGSC
    ngsc = shannon_entropy / max_entropy
    
    return ngsc

def analyze_choices():
    """
    Uses the NGSC calculation to analyze several answer choices.
    """
    print("--- Analyzing Answer Choices with Python ---")

    # --- Analysis for Choice M (NGSC=0) ---
    # Case: Maximum synchrony (one dominant component)
    # This corresponds to low complexity and high functional connectivity.
    eigs_max_sync = [100, 0.1, 0.05, 0.01]
    ngsc_max_sync = calculate_ngsc(eigs_max_sync)
    print(f"\n[Analysis for M]: If one component dominates (e.g., eigenvalues = {eigs_max_sync}), representing maximum synchrony:")
    print(f"The calculated NGSC is approximately {ngsc_max_sync:.4f}.")
    print("This shows that NGSC = 0 (or close to 0) means maximum SYNCHRONY, not independence. Choice M is incorrect.")

    # --- Analysis for Choice D (NGSC=1) ---
    # Case: Maximum desynchrony (all components are equal)
    # This corresponds to high complexity and low functional connectivity.
    eigs_max_desync = [25, 25, 25, 25]
    ngsc_max_desync = calculate_ngsc(eigs_max_desync)
    print(f"\n[Analysis for D]: If all components are equal (e.g., eigenvalues = {eigs_max_desync}), representing maximum desynchronization:")
    print(f"The calculated NGSC is {ngsc_max_desync:.4f}.")
    print("This shows that maximum complexity (NGSC=1) corresponds to LOW functional connectivity. Choice D is incorrect.")

    # --- Analysis for Choices E, H, N (Range of NGSC) ---
    # Case: A typical scenario with some distribution
    eigs_general = [70, 15, 10, 5]
    ngsc_general = calculate_ngsc(eigs_general)
    print(f"\n[Analysis for E, H, N]: For a general case (e.g., eigenvalues = {eigs_general}):")
    print(f"The calculated NGSC is {ngsc_general:.4f}.")
    print("The NGSC value is always between 0 and 1. Therefore, NGSC=2 (Choice E), NGSC<0 (Choice H), and NGSC=-1 (Choice N) are all mathematically impossible.")

if __name__ == '__main__':
    analyze_choices()