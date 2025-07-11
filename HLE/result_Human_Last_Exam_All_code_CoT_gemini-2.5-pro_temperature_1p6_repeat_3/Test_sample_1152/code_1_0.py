# Data from the experiments
mortality_non_infected = 10  # %
# Experiment 1 data for Fungus A
mortality_fungus_a = 35      # % (for most pollens)
# Experiment 3 data for Fungus B
mortality_fungus_b = 20      # %
# Experiment 4 data for Fungus C
mortality_fungus_c = 10      # %

print("Analyzing the nature of each fungus based on mortality rates...")
print("-" * 30)

# --- Analysis for Fungus A ---
is_pathogen_a = mortality_fungus_a > mortality_non_infected
print("Analysis for Fungus A:")
print(f"Is Fungus A mortality ({mortality_fungus_a}%) greater than non-infected mortality ({mortality_non_infected}%)?")
print(f"Result: {is_pathogen_a}")
if is_pathogen_a:
    print("Conclusion: Fungus A increases mortality, so it is a pathogen.\n")
else:
    print("Conclusion: Fungus A does not increase mortality.\n")


# --- Analysis for Fungus B ---
is_pathogen_b = mortality_fungus_b > mortality_non_infected
print("Analysis for Fungus B:")
print(f"Is Fungus B mortality ({mortality_fungus_b}%) greater than non-infected mortality ({mortality_non_infected}%)?")
print(f"Result: {is_pathogen_b}")
if is_pathogen_b:
    print("Conclusion: Fungus B increases mortality, so it is a pathogen.\n")
else:
    print("Conclusion: Fungus B does not increase mortality.\n")


# --- Analysis for Fungus C ---
# A commensal does not cause harm. The primary indicator of harm is increased mortality.
is_harmful_c = mortality_fungus_c > mortality_non_infected
print("Analysis for Fungus C:")
print(f"Does Fungus C increase mortality? Checking if {mortality_fungus_c}% > {mortality_non_infected}%.")
print(f"Result: {is_harmful_c}")
# We can also check if mortality is equal
is_commensal_c = mortality_fungus_c == mortality_non_infected
print(f"Is Fungus C mortality ({mortality_fungus_c}%) equal to non-infected mortality ({mortality_non_infected}%)?")
print(f"Result: {is_commensal_c}")
if is_commensal_c:
    print("Conclusion: Fungus C does not increase mortality. Furthermore, experiment 5 shows it often increases egg production. Therefore, it is not a pathogen and can be classified as a commensal (or mutualist).\n")
else:
    print("Conclusion: Fungus C affects mortality, so it is likely a pathogen.\n")

print("-" * 30)
print("Final Assessment:")
print("The evidence shows that Fungus A and B are pathogens, and Fungus C is a commensal.")
print("This directly corresponds to answer choice I.")
