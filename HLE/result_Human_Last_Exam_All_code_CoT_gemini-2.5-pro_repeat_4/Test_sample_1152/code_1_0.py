import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Data from the experiments ---
# Baseline mortality for non-infected bees
baseline_mortality = 10  # in percent

# Fungus A data
fungus_A_mortality = {
    "buck": 35, "sunflower": 10, "lavender": 35, "canola": 35,
    "milkweed": 35, "aster": 35, "mixed": 35
}
# Taking the common mortality rate for infected bees (excluding the protective sunflower pollen)
fungus_A_general_mortality = 35

# Fungus B data
fungus_B_mortality = {
    "buck": 20, "sunflower": 20, "lavender": 20, "canola": 20,
    "milkweed": 20, "aster": 20, "mixed": 20
}
# The mortality rate is uniform for all pollens
fungus_B_general_mortality = 20

# Fungus C data
fungus_C_mortality = {
    "buck": 10, "sunflower": 10, "lavender": 10, "canola": 10,
    "milkweed": 10, "aster": 10, "mixed": 10
}
# The mortality rate is uniform for all pollens
fungus_C_general_mortality = 10
fungus_C_productivity = {
    "buck": {"not_infected": 45, "infected": 60},
    "mixed": {"not_infected": 32, "infected": 52}
}

# --- Analysis ---

print("Step 1: Analyze if Fungus A and B are pathogens.")
print(f"The mortality rate of non-infected honeybees is {baseline_mortality}%.")

# Analysis of Fungus A
print("\nAnalysis of Fungus A:")
print(f"The mortality rate of honeybees infected with Fungus A (fed on most pollens) is {fungus_A_general_mortality}%.")
print(f"Comparing the mortality rates: {fungus_A_general_mortality}% (Fungus A) > {baseline_mortality}% (baseline).")
print("Since Fungus A increases the mortality rate, it is a pathogen.")

# Analysis of Fungus B
print("\nAnalysis of Fungus B:")
print(f"The mortality rate of honeybees infected with Fungus B is {fungus_B_general_mortality}%.")
print(f"Comparing the mortality rates: {fungus_B_general_mortality}% (Fungus B) > {baseline_mortality}% (baseline).")
print("Since Fungus B also increases the mortality rate, it is also a pathogen.")

print("\n----------------------------------------------------")

print("\nStep 2: Analyze if Fungus C is a commensal.")
print(f"The mortality rate of honeybees infected with Fungus C is {fungus_C_general_mortality}%.")
print(f"Comparing the mortality rates: {fungus_C_general_mortality}% (Fungus C) is equal to {baseline_mortality}% (baseline).")
print("Fungus C does not increase mortality. Let's check productivity.")

# Analysis of Fungus C Productivity
buck_prod = fungus_C_productivity["buck"]
print(f"For bees fed on buck pollen, egg production increased from {buck_prod['not_infected']} (not infected) to {buck_prod['infected']} (infected).")

mixed_prod = fungus_C_productivity["mixed"]
print(f"For bees fed on mixed pollen, egg production increased from {mixed_prod['not_infected']} (not infected) to {mixed_prod['infected']} (infected).")
print("Since Fungus C does not cause harm (no increased mortality) and can even have a beneficial effect (increased productivity), it is best described as a commensal.")

print("\n----------------------------------------------------")
print("\nConclusion:")
print("The analysis shows that Fungus A and B are pathogens, and Fungus C is a commensal.")
print("This corresponds to answer choice I.")


# --- Final Output ---
# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_str = captured_output.getvalue()
# Print the captured output
print(output_str)
print("<<<I>>>")