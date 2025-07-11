import collections

# Baseline mortality rate for non-infected honeybees
non_infected_mortality = 10

# Experiment 1 & 2 Data: Fungus A (Pathogen)
# Mortality rate with no effective treatment is 35%
mortality_fungus_A = 35

# Experiment 3 Data: Fungus B (Pathogen)
mortality_fungus_B = 20

# Experiment 4 & 5 Data: Fungus C (Commensal)
# Mortality rate is 10%, and productivity increases in many cases.
mortality_fungus_C = 10
eggs_C_uninfected_buck = 45
eggs_C_infected_buck = 60


print("Step 1: Determine if Fungus A is a pathogen.")
print(f"The mortality rate for bees infected with Fungus A is {mortality_fungus_A}%.")
print(f"The baseline mortality rate for non-infected bees is {non_infected_mortality}%.")
if mortality_fungus_A > non_infected_mortality:
    print(f"Since {mortality_fungus_A} > {non_infected_mortality}, Fungus A is a pathogen.\n")
else:
    print(f"Since {mortality_fungus_A} <= {non_infected_mortality}, Fungus A is not a pathogen.\n")


print("Step 2: Determine if Fungus B is a pathogen.")
print(f"The mortality rate for bees infected with Fungus B is {mortality_fungus_B}%.")
if mortality_fungus_B > non_infected_mortality:
    print(f"Since {mortality_fungus_B} > {non_infected_mortality}, Fungus B is a pathogen.\n")
else:
    print(f"Since {mortality_fungus_B} <= {non_infected_mortality}, Fungus B is not a pathogen.\n")


print("Step 3: Determine if Fungus C is a pathogen or commensal.")
print(f"The mortality rate for bees infected with Fungus C is {mortality_fungus_C}%.")
if mortality_fungus_C > non_infected_mortality:
    print(f"Since {mortality_fungus_C} > {non_infected_mortality}, Fungus C is a pathogen.")
else:
    print(f"Since {mortality_fungus_C} is not greater than {non_infected_mortality}, Fungus C does not increase mortality.")
    print("Let's check its effect on productivity:")
    print(f"For Buck pollen, egg count increased from {eggs_C_uninfected_buck} (uninfected) to {eggs_C_infected_buck} (infected).")
    print("Since it does not cause harm and can be beneficial, it is not a pathogen. 'Commensal' is a suitable term.\n")

print("Conclusion: The evidence supports that Fungus A and B are pathogens, and Fungus C is a commensal.")
print("This matches answer choice I.")
print("<<<I>>>")