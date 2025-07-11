import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# --- Data Representation ---
baseline_mortality = 10

# Experiment 1 & 2 Data: Fungus A
fungus_A_data = {
    "buck": {"mortality": 35, "eggs_i": 10},
    "sunflower": {"mortality": 10, "eggs_i": 20},
}

# Experiment 3 Data: Fungus B
fungus_B_data = {
    "buck": {"mortality": 20},
    "sunflower": {"mortality": 20},
}

# Experiment 4 & 5 Data: Fungus C
fungus_C_data = {
    "buck": {"mortality": 10, "eggs_ni": 45, "eggs_i": 60},
    "sunflower": {"mortality": 10, "eggs_ni": 30, "eggs_i": 25},
}


# --- Step-by-step Analysis ---
print("Analysis of Experimental Data:\n")

# 1. Analyze Fungus A
mortality_A = fungus_A_data["buck"]["mortality"]
is_A_pathogen = mortality_A > baseline_mortality
print(f"Step 1: Analyze Fungus A's nature.")
print(f"  - The highest mortality rate with Fungus A infection is {mortality_A}%.")
print(f"  - The baseline mortality for non-infected bees is {baseline_mortality}%.")
print(f"  - Since {mortality_A} > {baseline_mortality}, Fungus A increases mortality and is a pathogen.")
print("-" * 20)

# 2. Analyze Fungus B
mortality_B = fungus_B_data["buck"]["mortality"]
is_B_pathogen = mortality_B > baseline_mortality
print(f"Step 2: Analyze Fungus B's nature.")
print(f"  - The mortality rate with Fungus B infection is {mortality_B}%.")
print(f"  - The baseline mortality is {baseline_mortality}%.")
print(f"  - Since {mortality_B} > {baseline_mortality}, Fungus B increases mortality and is a pathogen.")
print("-" * 20)

# 3. Analyze Fungus C
mortality_C = fungus_C_data["buck"]["mortality"]
is_C_commensal = mortality_C == baseline_mortality
print(f"Step 3: Analyze Fungus C's nature.")
print(f"  - The mortality rate with Fungus C infection is {mortality_C}%.")
print(f"  - This is equal to the baseline mortality of {baseline_mortality}%.")
print(f"  - Since it does not cause additional harm (death), Fungus C is not a pathogen. It is a commensal.")
print(f"  - In fact, with buck pollen, egg production increases from {fungus_C_data['buck']['eggs_ni']} to {fungus_C_data['buck']['eggs_i']}, suggesting a beneficial relationship.")
print("-" * 20)

# 4. Final Conclusion based on analysis
print("Conclusion:")
print("Based on the analysis of mortality rates:")
print(f"  - Fungus A is a pathogen: {is_A_pathogen}")
print(f"  - Fungus B is a pathogen: {is_B_pathogen}")
print(f"  - Fungus C is a commensal: {is_C_commensal}")
print("\nThe statement 'Fungus A and B are pathogens. Fungus C is a commensal.' is fully supported by the data.")
print("This corresponds to answer choice I.")


# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output_string = captured_output.getvalue()

# Print the final captured output
print(output_string)
print("<<<I>>>")