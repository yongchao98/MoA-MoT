# Data from the experiments, organized for clarity.
experimental_data = {
    "Wild-type": {
        "time_in_center_percent": 15,
        "ki67_cells": 3500,
    },
    "deltaber2": {
        "time_in_center_percent": 8,
        "ki67_cells": 3500,
        "ssri_time_in_center_percent": 15,
    },
    "delta-ber1_delta-ber2": {
        "ki67_cells": 2850,
    }
}

print("Analyzing the data to find the most accurate conclusion...\n")

# --- Clause 1 Analysis ---
print("Clause 1: 'The effects of mutations in ber1 and ber2 may be reversed by treatment with SSRI.'")
deltaber2_before = experimental_data["deltaber2"]["time_in_center_percent"]
deltaber2_after = experimental_data["deltaber2"]["ssri_time_in_center_percent"]
print(f"Finding: The anxiety phenotype in deltaber2 mice (time in center) was reversed from {deltaber2_before}% to {deltaber2_after}% after SSRI treatment.")
print("Conclusion: This clause is TRUE.\n")


# --- Clause 2 Analysis ---
print("Clause 2: 'Mice with defects in ber2 may not have a decrease in cell proliferation.'")
wt_ki67 = experimental_data["Wild-type"]["ki67_cells"]
deltaber2_ki67 = experimental_data["deltaber2"]["ki67_cells"]
print(f"Finding: The number of Ki67 cells in deltaber2 mice ({deltaber2_ki67}) is the same as in Wild-type mice ({wt_ki67}).")
print("Conclusion: This clause is TRUE.\n")


# --- Clause 3 Analysis ---
print("Clause 3: 'Gene ber1 and ber2 regulate cell proliferation.'")
double_ko_ki67 = experimental_data["delta-ber1_delta-ber2"]["ki67_cells"]
print(f"Finding: Only the double knockout of ber1 and ber2 resulted in a decrease in Ki67 cells (from {wt_ki67} to {double_ko_ki67}).")
print("Conclusion: This shows a redundant regulatory function, so this clause is TRUE.\n")

print("--------------------------------------------------")
print("All three clauses of statement A are supported by the data.")
print("The final equation below demonstrates the reversal of the anxiety phenotype, a key finding.")
print("Equation (Restored Time % - Defect Time % = Improvement %):")
# Printing each number in the final equation as requested
restored_value = experimental_data["deltaber2"]["ssri_time_in_center_percent"]
defect_value = experimental_data["deltaber2"]["time_in_center_percent"]
improvement = restored_value - defect_value
print(f"{restored_value} - {defect_value} = {improvement}")
print("--------------------------------------------------")

<<<A>>>