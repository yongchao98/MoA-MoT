# Baseline data
baseline_mortality_rate = 10

# Experiment 1 & 4 data for Fungus A and C
fungus_A_max_mortality = 35
fungus_C_mortality = 10

# Experiment 3 data for Fungus B
fungus_B_mortality = 20

# Experiment 2 & 5 productivity data (infected bees)
# Format: { 'Pollen': eggs }
productivity_infected_A = {
    'buck': 10,
    'sunflower': 20
}

productivity_infected_C = {
    'sunflower': 25,
    'lavender': 50
}

# --- Analysis ---

print("Analyzing the classification of each fungus based on the data provided.")
print("-" * 60)

# Analyze Fungus A
is_A_pathogen = fungus_A_max_mortality > baseline_mortality_rate
print(f"1. Fungus A has a maximum mortality rate of {fungus_A_max_mortality}%.")
print(f"   The baseline mortality for non-infected bees is {baseline_mortality_rate}%.")
if is_A_pathogen:
    print(f"   Since {fungus_A_max_mortality} > {baseline_mortality_rate}, Fungus A is a pathogen.")
else:
    print(f"   Since {fungus_A_max_mortality} <= {baseline_mortality_rate}, Fungus A is not a pathogen.")

print("-" * 60)

# Analyze Fungus B
is_B_pathogen = fungus_B_mortality > baseline_mortality_rate
print(f"2. Fungus B has a mortality rate of {fungus_B_mortality}%.")
print(f"   The baseline mortality for non-infected bees is {baseline_mortality_rate}%.")
if is_B_pathogen:
    print(f"   Since {fungus_B_mortality} > {baseline_mortality_rate}, Fungus B is a pathogen.")
else:
    print(f"   Since {fungus_B_mortality} <= {baseline_mortality_rate}, Fungus B is not a pathogen.")

print("-" * 60)

# Analyze Fungus C
is_C_pathogen = fungus_C_mortality > baseline_mortality_rate
print(f"3. Fungus C has a mortality rate of {fungus_C_mortality}%.")
print(f"   The baseline mortality for non-infected bees is {baseline_mortality_rate}%.")
if is_C_pathogen:
    print(f"   Since {fungus_C_mortality} > {baseline_mortality_rate}, Fungus C is a pathogen.")
else:
    # A commensal does not harm the host. Increased productivity is a benefit.
    print(f"   Since {fungus_C_mortality} is not greater than {baseline_mortality_rate} and it increases egg production, Fungus C is a commensal.")

print("-" * 60)
print("\nFinal Conclusion based on the analysis:")

if is_A_pathogen and is_B_pathogen and not is_C_pathogen:
    print("Fungus A and B are pathogens. Fungus C is a commensal.")
else:
    print("The analysis does not match the proposed conclusion.")

<<<I>>>