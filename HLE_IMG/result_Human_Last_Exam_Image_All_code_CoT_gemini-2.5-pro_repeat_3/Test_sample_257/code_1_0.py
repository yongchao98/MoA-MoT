def get_splitting_pattern_name(num_peaks):
    """Converts the number of peaks into the common name for the multiplet."""
    patterns = {
        1: "singlet",
        2: "doublet",
        3: "triplet",
        4: "quartet",
        5: "quintet",
        6: "sextet",
    }
    return patterns.get(num_peaks, f"{num_peaks}-peak multiplet")

# Based on chemical analysis, we identify the properties of the most deshielded proton.

# For integration: determine the number of protons giving rise to the signal.
# The most deshielded proton is the unique central H on the top ring. There is only one.
num_protons_in_signal = 1

# For splitting: determine the number of equivalent neighboring protons ('n').
# This proton has two neighbors on the adjacent carbons. Due to symmetry, they are equivalent.
n_neighbors = 2

# Apply the n+1 rule to find the number of peaks in the signal.
num_peaks = n_neighbors + 1

# Get the name of the splitting pattern from the number of peaks.
splitting_pattern = get_splitting_pattern_name(num_peaks)

# Print the detailed analysis and the final answer.
print("Analysis of the 1H NMR signal for the most deshielded proton:")
print("-" * 50)

print("1. Integration:")
print(f"The signal corresponds to a single, unique proton in the molecule.")
print(f"   Integration = {num_protons_in_signal}H")
print("-" * 50)

print("2. Splitting Pattern (Multiplicity):")
print("The multiplicity is determined by the n+1 rule.")
print(f"   Number of equivalent neighboring protons (n) = {n_neighbors}")
print(f"   Number of peaks = n + 1 = {n_neighbors} + 1 = {num_peaks}")
print(f"   A signal with {num_peaks} peaks is called a {splitting_pattern}.")
print("-" * 50)

print(f"Conclusion: The highest deshielded proton peak is a {splitting_pattern} with an integration of {num_protons_in_signal}H.")