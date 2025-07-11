import numpy as np

def explain_chemotype_generation(group_name, n_specimens, n_peaks):
    """
    This function demonstrates the logic from the R script to explain
    why only one chemotype is generated per group.
    """
    print(f"--- Analyzing the '{group_name}' group generation ---")
    
    # In the R script, a single baseline vector is created for the entire group.
    # Let's simulate that.
    baseline = np.random.uniform(0, 1, n_peaks)
    print(f"1. A single, common baseline vector of size {n_peaks} is created for all {n_specimens} specimens.")

    # The script then subtracts this same baseline from every specimen in the group.
    print("2. This exact same baseline is subtracted from every specimen in the group.")
    
    # The final step is normalization. This rescales the data but does not create new clusters.
    print("3. Each specimen is then normalized.")

    # The conclusion is drawn from this process.
    print(f"\nConclusion for '{group_name}': Because all specimens are processed identically with a common baseline,")
    print("they are not separated into distinct subgroups. They will form a single data cluster.")
    
    # The number of chemotypes is the number of clusters.
    num_chemotypes = 1
    print(f"This corresponds to {num_chemotypes} chemotype for the {group_name} group.")
    return num_chemotypes

# --- Main Execution ---
print("Based on the logic of the provided R script, we can determine the expected number of chemotypes.\n")

# Determine chemotypes for the control group
control_chemotypes = explain_chemotype_generation("control", 360, 340)

print("\n" + "="*60 + "\n")

# Determine chemotypes for the water stress group
water_stress_chemotypes = explain_chemotype_generation("water_stress", 510, 340)

print("\n" + "="*60 + "\n")
print("Final Answer:")
print(f"The number of expected chemotypes for the control group is: {control_chemotypes}")
print(f"The number of expected chemotypes for the water stressed group is: {water_stress_chemotypes}")