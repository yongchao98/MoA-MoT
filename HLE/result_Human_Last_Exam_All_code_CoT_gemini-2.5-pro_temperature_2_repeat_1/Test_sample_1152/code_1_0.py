# Data from the experiments
# Mortality rates (%)
mortality_non_infected = 10
# From Experiment 1: Most pollen types showed a 35% mortality rate with Fungus A
mortality_fungus_A = 35
# From Experiment 3: All pollen types showed a 20% mortality rate with Fungus B
mortality_fungus_B = 20
# From Experiment 4: All pollen types showed a 10% mortality rate with Fungus C
mortality_fungus_C = 10

# Productivity (number of eggs) for a representative pollen (e.g., Lavender)
# From Experiment 2 and 5
productivity_lavender_not_infected = 30
# From Experiment 5
productivity_lavender_infected_C = 50

print("--- Analysis of Fungi ---")

# Step 1: Analyze Fungus A
print(f"Mortality of non-infected bees: {mortality_non_infected}%")
print(f"Mortality of bees infected with Fungus A: {mortality_fungus_A}%")
if mortality_fungus_A > mortality_non_infected:
    print("Conclusion: Fungus A causes increased mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus A does not increase mortality.")

print("-" * 20)

# Step 2: Analyze Fungus B
print(f"Mortality of non-infected bees: {mortality_non_infected}%")
print(f"Mortality of bees infected with Fungus B: {mortality_fungus_B}%")
if mortality_fungus_B > mortality_non_infected:
    print("Conclusion: Fungus B causes increased mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus B does not increase mortality.")

print("-" * 20)

# Step 3: Analyze Fungus C
print(f"Mortality of non-infected bees: {mortality_non_infected}%")
print(f"Mortality of bees infected with Fungus C: {mortality_fungus_C}%")
if mortality_fungus_C > mortality_non_infected:
    print("Conclusion: Fungus C increases mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus C does not increase mortality.")
    # Check productivity impact
    print(f"\nProductivity of non-infected bees (Lavender pollen): {productivity_lavender_not_infected} eggs.")
    print(f"Productivity of bees infected with Fungus C (Lavender pollen): {productivity_lavender_infected_C} eggs.")
    if productivity_lavender_infected_C > productivity_lavender_not_infected:
        print("Productivity increased with Fungus C infection. This suggests Fungus C is not harmful and may be beneficial.")
    print("Conclusion: Since Fungus C does not cause harm and can increase productivity, it is best described as a commensal.")

print("\n--- Final Evaluation ---")
print("Based on the analysis:")
print("- Fungus A is a pathogen.")
print("- Fungus B is a pathogen.")
print("- Fungus C is a commensal.")
print("The statement that accurately reflects this is: 'Fungus A and B are pathogens. Fungus C is a commensal.'")

<<<I>>>