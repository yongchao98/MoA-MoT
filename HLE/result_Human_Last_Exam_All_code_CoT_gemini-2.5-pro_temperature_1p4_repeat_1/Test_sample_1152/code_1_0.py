import sys
# Redirect stdout to a variable to control output format
# This is a trick to make sure we only have one code block.
original_stdout = sys.stdout
from io import StringIO
sys.stdout = captured_output = StringIO()


# --- Data from the experiments ---
baseline_mortality = 10  # Mortality rate of non-infected bees in %

# Experiment 1 & 3 & 4: Mortality data for different fungi
# We use the general mortality rate observed for the fungus
fungus_A_mortality = 35
fungus_B_mortality = 20
fungus_C_mortality = 10

# Experiment 5: Productivity data for Fungus C infection
# (number of eggs)
productivity_not_infected = {
    'buck': 45,
    'lavender': 30,
    'canola': 30
}
productivity_C_infected = {
    'buck': 60,
    'lavender': 50,
    'canola': 50
}


# --- Analysis ---
print("Step-by-step analysis to determine the nature of each fungus:\n")

# Analyze Fungus A
print("--- Analysis of Fungus A ---")
print(f"The baseline mortality rate of honeybees is {baseline_mortality}%.")
print(f"The mortality rate for bees infected with Fungus A is {fungus_A_mortality}%.")
if fungus_A_mortality > baseline_mortality:
    print(f"Comparison: {fungus_A_mortality} > {baseline_mortality}")
    print("Conclusion: Fungus A increases mortality, so it is a pathogen.\n")
else:
    print(f"Comparison: {fungus_A_mortality} <= {baseline_mortality}")
    print("Conclusion: Fungus A does not increase mortality.\n")


# Analyze Fungus B
print("--- Analysis of Fungus B ---")
print(f"The baseline mortality rate of honeybees is {baseline_mortality}%.")
print(f"The mortality rate for bees infected with Fungus B is {fungus_B_mortality}%.")
if fungus_B_mortality > baseline_mortality:
    print(f"Comparison: {fungus_B_mortality} > {baseline_mortality}")
    print("Conclusion: Fungus B increases mortality, so it is a pathogen.\n")
else:
    print(f"Comparison: {fungus_B_mortality} <= {baseline_mortality}")
    print("Conclusion: Fungus B does not increase mortality.\n")


# Analyze Fungus C
print("--- Analysis of Fungus C ---")
# Part 1: Mortality
print("Part 1: Analyzing mortality rate.")
print(f"The baseline mortality rate of honeybees is {baseline_mortality}%.")
print(f"The mortality rate for bees infected with Fungus C is {fungus_C_mortality}%.")
if fungus_C_mortality > baseline_mortality:
    print(f"Comparison: {fungus_C_mortality} > {baseline_mortality}")
    print("Intermediate Conclusion: Fungus C appears to be a pathogen.\n")
    is_pathogen = True
else:
    print(f"Comparison: {fungus_C_mortality} <= {baseline_mortality}")
    print("Intermediate Conclusion: Fungus C does not increase mortality.\n")
    is_pathogen = False

# Part 2: Productivity if not a pathogen
if not is_pathogen:
    print("Part 2: Analyzing productivity since mortality did not increase.")
    pollen_type = 'buck'
    eggs_not_infected = productivity_not_infected[pollen_type]
    eggs_infected = productivity_C_infected[pollen_type]
    print(f"For {pollen_type} pollen, non-infected bees produced {eggs_not_infected} eggs.")
    print(f"For {pollen_type} pollen, bees infected with Fungus C produced {eggs_infected} eggs.")
    if eggs_infected > eggs_not_infected:
        print(f"Comparison: {eggs_infected} > {eggs_not_infected}")
        print("Conclusion: Fungus C does not cause harm and increases productivity, so it is a commensal.\n")
    else:
        print(f"Comparison: {eggs_infected} <= {eggs_not_infected}")
        print("Conclusion: Fungus C does not increase productivity in this case.\n")

# --- Final Summary ---
print("====================================")
print("Final Summary based on the analysis:")
print("- Fungus A is a pathogen.")
print("- Fungus B is a pathogen.")
print("- Fungus C is a commensal.")
print("This corresponds to answer choice I.")
print("====================================")

# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())
<<<I>>>