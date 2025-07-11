import numpy as np

def calculate_ngsc(eigenvalues):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC).
    NGSC is the normalized entropy of normalized eigenvalues.
    The function demonstrates that the output is always in the range [0, 1].
    """
    # Eigenvalues represent variance and must be non-negative.
    eigenvalues = np.array(eigenvalues, dtype=float)
    if np.any(eigenvalues < 0):
        raise ValueError("Eigenvalues cannot be negative.")
        
    total_variance = np.sum(eigenvalues)
    if total_variance == 0:
        # If there is no variance, complexity is undefined or zero.
        return 0

    # Normalize eigenvalues to unit sum to create a probability distribution p_i.
    normalized_eigenvalues = eigenvalues / total_variance
    
    # Filter out zero values to avoid log(0) which is -inf.
    # p_i * log(p_i) for p_i=0 is 0.
    non_zero_eigs = normalized_eigenvalues[normalized_eigenvalues > 0]
    
    # The number of components, N, is the total number of eigenvalues.
    N = len(eigenvalues)
    if N <= 1:
        # For a single component, entropy is 0.
        return 0
        
    # Calculate Shannon entropy H(p) = -sum(p_i * log(p_i)).
    # We use the natural logarithm (base e).
    # The sum is over non_zero_eigs, as 0 * log(0) is 0.
    entropy = -np.sum(non_zero_eigs * np.log(non_zero_eigs))
    
    # The maximum possible entropy for a distribution with N outcomes is log(N).
    max_entropy = np.log(N)
    if max_entropy == 0:
        return 0
    
    # Normalize the calculated entropy by the maximum possible entropy.
    ngsc = entropy / max_entropy
    
    return ngsc

# --- Analysis of Answer Choice E ---
# E. If NGSC = 2, then the Gini coefficient ... will be 0.5
print("--- Analyzing Answer Choice E with Python ---")
print("Statement: 'If NGSC = 2, then the Gini coefficient ... will be 0.5'")
print("\nFirst, we test the premise 'If NGSC = 2'. We will check the range of NGSC with two extreme theoretical cases.")

# The number of components (e.g., voxels or parcels) is N. Let's use 50.
N = 50
print(f"Let's assume a system with N = {N} principal components.\n")

# Case 1: Minimum Complexity (perfect synchrony).
# All variance is explained by one component. Eigenvalues are [1, 0, 0, ...].
min_complexity_eigs = np.zeros(N)
min_complexity_eigs[0] = 1.0
ngsc_min = calculate_ngsc(min_complexity_eigs)

print(f"Case 1: Minimum Complexity (Perfect Synchrony)")
print(f"Eigenvalue distribution is [1.0, 0.0, 0.0, ...]")
print(f"The calculated NGSC is: {ngsc_min:.4f}")
print("This is the absolute minimum NGSC, representing zero entropy.\n")


# Case 2: Maximum Complexity (perfect desynchrony).
# All variance is distributed equally among all components. Eigenvalues are [1/N, 1/N, ...].
max_complexity_eigs = np.ones(N) * (1.0 / N)
ngsc_max = calculate_ngsc(max_complexity_eigs)

print(f"Case 2: Maximum Complexity (Perfect Desynchrony)")
print(f"Eigenvalue distribution is [1/{N}, 1/{N}, ...]")
print(f"The calculated NGSC is: {ngsc_max:.4f}")
print("This is the absolute maximum NGSC, representing maximum entropy.\n")

# Conclusion from the code
print("--- Conclusion ---")
print("As demonstrated, NGSC is a normalized value, mathematically bound to the range [0, 1].")
print("An NGSC value of 2 is therefore impossible.")
print("Because the premise of statement E is false, the entire statement is logically incorrect.")
