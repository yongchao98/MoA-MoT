import numpy as np

def generate_chemistry_py(n_specimens, n_peaks):
    """
    A Python implementation of the R script's core logic.
    This function simulates one group of specimens.
    """
    # Generate initial peak data with random variation for each specimen.
    # In R: lapply(1:n_peaks,\(m) rnorm(n_specimens, m, 1)) |> do.call(cbind, args = _)
    peaks = np.zeros((n_specimens, n_peaks))
    for m in range(n_peaks):
        peaks[:, m] = np.random.normal(loc=m + 1, scale=1, size=n_specimens)

    # A single baseline is generated for the entire group.
    # In R: baseline = runif(n_peaks, 0, 1)
    baseline = np.random.uniform(0, 1, n_peaks)

    # This single baseline is subtracted from all specimens in the group.
    # This means all specimens are systematically adjusted in the exact same way.
    # In R: apply(peaks[t==treatment,], 1, \(x)x-baseline)
    peaks_adjusted = peaks - baseline
    
    # The R script removes the first specimen for an unknown reason. We replicate that.
    # In R: return(peaks[-1,])
    final_peaks = peaks_adjusted[1:, :]

    return final_peaks

def explain_chemotypes():
    """
    Explains the expected number of chemotypes based on the data generation logic.
    """
    print("Analyzing the data generation process:")
    print("-" * 40)
    
    # Analyze Control Group
    print("1. For the 'control' group:")
    print("   - A single function call generates all 360 specimens.")
    print("   - A single 'baseline' vector is created.")
    print("   - This SAME baseline is subtracted from ALL 360 control specimens.")
    print("   - Conclusion: All control specimens belong to one systematically-defined group.")
    control_chemotypes = 1
    print(f"   - Expected number of chemotypes for the control group = {control_chemotypes}\n")

    # Analyze Water Stress Group
    print("2. For the 'water_stress' group:")
    print("   - A SEPARATE function call generates all 510 specimens.")
    print("   - A NEW, single 'baseline' vector is created for this group.")
    print("   - This SAME new baseline is subtracted from ALL 510 water-stressed specimens.")
    print("   - Conclusion: All water-stressed specimens belong to one systematically-defined group.")
    water_stress_chemotypes = 1
    print(f"   - Expected number of chemotypes for the water_stress group = {water_stress_chemotypes}\n")
    
    print("-" * 40)
    print("Final Answer based on the script's logic:")
    print(f"Control Group Chemotypes: {control_chemotypes}")
    print(f"Water Stress Group Chemotypes: {water_stress_chemotypes}")

# Running the explanation
# NOTE: The generate_chemistry_py function is provided to show the logic, 
# but we call explain_chemotypes() to get the direct answer.
explain_chemotypes()
