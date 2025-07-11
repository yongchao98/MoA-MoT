import numpy as np

def generate_chemistry(n_specimens, n_peaks):
    """
    A Python implementation of the R script's data generation logic.
    For a single call, it generates a single population of samples.
    """
    # Create peaks where each peak (column) has a unique mean.
    # This corresponds to R's `lapply(1:n_peaks, \(m) rnorm(n_specimens, m, 1))`
    means = np.arange(1, n_peaks + 1)
    peaks = np.random.normal(loc=means, scale=1, size=(n_specimens, n_peaks))

    # Generate a single baseline vector for this entire batch.
    # This corresponds to R's `baseline = runif(n_peaks, 0, 1)`
    baseline = np.random.uniform(0, 1, n_peaks)

    # In the R script, because it is called with a single treatment name,
    # the same baseline is subtracted from all specimens in the batch.
    # This is a simple broadcasted subtraction in NumPy.
    peaks_baselined = peaks - baseline
    
    # The R script oddly removes the first row. We will ignore this minor
    # detail as it does not change the number of populations (chemotypes).
    # return peaks_baselined[1:,:]
    return peaks_baselined

# --- Main Analysis ---

# In the R script, the control group is generated in one self-contained call.
# All generated specimens come from one set of parameters (one baseline).
# Therefore, they constitute a single population or chemotype.
num_chemotypes_control = 1

# The water stress group is generated in a second, separate call.
# All specimens from this call also come from one set of parameters.
# Therefore, they also constitute a single population or chemotype.
num_chemotypes_stress = 1

print("Analyzing the data generation process from the script:")
print(f"The 'control' group is created from a single, homogeneous process.")
print(f"This results in an expected number of chemotypes for the control group: {num_chemotypes_control}")
print("")
print(f"The 'water_stress' group is also created from a single, homogeneous process.")
print(f"This results in an expected number of chemotypes for the water stress group: {num_chemotypes_stress}")
print("")
print(f"Conclusion: We expect {num_chemotypes_control} chemotype for controls and {num_chemotypes_stress} for water-stressed.")
