# Data from the experiments
# Mortality Rates (%)
mortality_non_infected = 10
mortality_fungus_A = 35  # For most pollens
mortality_fungus_B = 20
mortality_fungus_C = 10

# Productivity data (eggs) for Fungus C experiment
productivity_buck_not_infected = 45
productivity_buck_infected_C = 60

# Step 1 & 2: Evaluate if Fungi are pathogens by comparing mortality rates
print("Step 1: Evaluating if Fungus A, B, and C are pathogens.")
print(f"The mortality rate of non-infected bees is {mortality_non_infected}%.")

# Analysis for Fungus A
print("\n--- Analysis of Fungus A ---")
print(f"The mortality rate of bees infected with Fungus A is {mortality_fungus_A}%.")
is_A_pathogen = mortality_fungus_A > mortality_non_infected
print(f"Is the mortality rate for Fungus A ({mortality_fungus_A}%) > the non-infected rate ({mortality_non_infected}%)? {is_A_pathogen}.")
if is_A_pathogen:
    print("Conclusion: Fungus A is a pathogen because it increases the mortality rate.")

# Analysis for Fungus B
print("\n--- Analysis of Fungus B ---")
print(f"The mortality rate of bees infected with Fungus B is {mortality_fungus_B}%.")
is_B_pathogen = mortality_fungus_B > mortality_non_infected
print(f"Is the mortality rate for Fungus B ({mortality_fungus_B}%) > the non-infected rate ({mortality_non_infected}%)? {is_B_pathogen}.")
if is_B_pathogen:
    print("Conclusion: Fungus B is a pathogen because it increases the mortality rate.")

# Analysis for Fungus C
print("\n--- Analysis of Fungus C ---")
print(f"The mortality rate of bees infected with Fungus C is {mortality_fungus_C}%.")
is_C_pathogen = mortality_fungus_C > mortality_non_infected
print(f"Is the mortality rate for Fungus C ({mortality_fungus_C}%) > the non-infected rate ({mortality_non_infected}%)? {is_C_pathogen}.")
if not is_C_pathogen:
    print("Conclusion: Fungus C does not increase the mortality rate, so it is not a pathogen.")

# Step 3: Classify Fungus C
print("\nStep 2: Classifying the relationship of Fungus C with honeybees.")
print("An organism that lives with a host without causing harm is a commensal.")
print(f"Since Fungus C infection results in a mortality rate of {mortality_fungus_C}%, which is equal to the non-infected rate, it does not harm the bees' survival.")
# Check productivity impact
productivity_change_buck = productivity_buck_infected_C - productivity_buck_not_infected
print(f"Furthermore, for bees fed buck pollen, Fungus C infection increased egg production from {productivity_buck_not_infected} to {productivity_buck_infected_C}, an increase of {productivity_change_buck} eggs.")
print("This lack of harm and potential benefit supports classifying Fungus C as a commensal (or even mutualist).")

# Step 4: Final Conclusion
print("\n--- Final Conclusion ---")
print("Based on the analysis:")
print("- Fungus A and B are pathogens.")
print("- Fungus C is a commensal.")
print("The statement that accurately reflects this is 'I. Fungus A and B are pathogens. Fungus C is a commensal.'")
<<<I>>>