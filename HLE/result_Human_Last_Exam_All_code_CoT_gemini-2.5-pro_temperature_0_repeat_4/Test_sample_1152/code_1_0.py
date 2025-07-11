# Step 1: Define the mortality rates from the experiments.
# The mortality rate of non-infected honeybees serves as our baseline.
baseline_mortality = 10  # in percent

# Mortality rates for bees infected with each fungus.
# We use the general rate, as some pollens can have a protective effect.
mortality_A = 35  # Mortality rate with Fungus A (e.g., on buck pollen)
mortality_B = 20  # Mortality rate with Fungus B
mortality_C = 10  # Mortality rate with Fungus C

# Step 2: Calculate the increase in mortality caused by each fungus.
increase_A = mortality_A - baseline_mortality
increase_B = mortality_B - baseline_mortality
increase_C = mortality_C - baseline_mortality

# Step 3: Analyze and classify each fungus.
print("Analysis of Fungi Based on Mortality Increase:")
print("="*50)

# Analysis for Fungus A
print("Fungus A Analysis:")
print(f"Mortality increase = Mortality with Fungus A - Baseline Mortality")
print(f"Mortality increase = {mortality_A}% - {baseline_mortality}% = {increase_A}%")
if increase_A > 0:
    print("Conclusion: Fungus A increases mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus A does not increase mortality, so it is not a pathogen.")
print("-" * 50)

# Analysis for Fungus B
print("Fungus B Analysis:")
print(f"Mortality increase = Mortality with Fungus B - Baseline Mortality")
print(f"Mortality increase = {mortality_B}% - {baseline_mortality}% = {increase_B}%")
if increase_B > 0:
    print("Conclusion: Fungus B increases mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus B does not increase mortality, so it is not a pathogen.")
print("-" * 50)

# Analysis for Fungus C
print("Fungus C Analysis:")
print(f"Mortality increase = Mortality with Fungus C - Baseline Mortality")
print(f"Mortality increase = {mortality_C}% - {baseline_mortality}% = {increase_C}%")
if increase_C > 0:
    print("Conclusion: Fungus C increases mortality, so it is a pathogen.")
else:
    print("Conclusion: Fungus C does not increase mortality, so it is a commensal.")
    # Experiment 5 shows that with some pollens, productivity even increases with Fungus C,
    # further supporting that it is not harmful.
print("="*50)

print("\nFinal Summary: The analysis shows that Fungus A and Fungus B are pathogens, while Fungus C is a commensal.")
<<<I>>>