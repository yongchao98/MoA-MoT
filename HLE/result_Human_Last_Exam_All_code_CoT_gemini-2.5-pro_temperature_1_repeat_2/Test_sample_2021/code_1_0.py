# Step 1: Define variables based on the problem description.
# Genetic parameters
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
# For a self-pollination cross (Ee x Ee), the offspring are 1/4 EE, 1/2 Ee, and 1/4 ee.
# The resistant phenotypes are Ee and ee, so the fraction is 1/2 + 1/4 = 3/4.
resistant_fraction_from_selfing = 3/4
# For a cross-pollination with the wild-type population (Ee x EE), the offspring are 1/2 EE and 1/2 Ee.
# The resistant phenotype is Ee, so the fraction is 1/2.
resistant_fraction_from_crossing = 1/2

# Molecular biology parameters
insertion_length_bp = 105
nucleotides_per_amino_acid = 3
avg_amino_acid_mass_daltons = 110
daltons_to_kda_conversion = 1000

# Step 2: Perform calculations and print the logical steps.

print("--- Analysis Part 1: Calculating Percentage of Resistant Offspring ---")
# The total fraction of resistant offspring is the sum of contributions from self-pollination and cross-pollination.
contribution_from_selfing = self_pollination_rate * resistant_fraction_from_selfing
contribution_from_crossing = cross_pollination_rate * resistant_fraction_from_crossing
total_resistant_fraction = contribution_from_selfing + contribution_from_crossing
total_resistant_percentage = total_resistant_fraction * 100

print("The calculation for the total fraction of resistant offspring is:")
print(f"({self_pollination_rate} * {resistant_fraction_from_selfing}) + ({cross_pollination_rate} * {resistant_fraction_from_crossing}) = {total_resistant_fraction}")
print(f"This corresponds to a theoretical percentage of {total_resistant_percentage:.2f}% resistant offspring.")
print("-" * 50)

print("--- Analysis Part 2: Determining Protein Activity and Interaction ---")
print("1. Protein Activity Analysis:")
print("   - With E3ub-wt, the level of Par22 decreased (700 -> 200 units).")
print("   - Conclusion: E3ub-wt is an active E3 ubiquitin ligase that targets Par22.")
print("   - With E3ub-insert105, the level of Par22 was not reduced (700 -> 3000 units).")
print("   - Conclusion: E3ub-insert105 is NOT an active E3 ubiquitin ligase.")
print("\n2. Protein Interaction Analysis:")
print("   - Par22 (50kDa) + E3ub-wt (60kDa) produced a single peak at 110kDa.")
print("   - Conclusion: E3ub-wt interacts with Par22 to form a complex.")
print("   - Par22 + E3ub-insert105 produced separate peaks for each protein.")
print("   - Conclusion: E3ub-insert105 does NOT interact with Par22.")
print("-" * 50)


print("--- Analysis Part 3: Calculating Protein Mass Increase ---")
# The mass increase is calculated from the number of amino acids added by the insertion.
amino_acids_added = insertion_length_bp / nucleotides_per_amino_acid
mass_increase_daltons = amino_acids_added * avg_amino_acid_mass_daltons
mass_increase_kda = mass_increase_daltons / daltons_to_kda_conversion

print(f"The {insertion_length_bp} nucleotide insertion adds {int(amino_acids_added)} amino acids to the protein.")
print("The calculation for the approximate mass increase is:")
print(f"({insertion_length_bp} nucleotides / {nucleotides_per_amino_acid}) * {avg_amino_acid_mass_daltons} Da = {mass_increase_daltons} Da")
print(f"This is approximately {mass_increase_kda:.1f} kDa.")
print("-" * 50)

print("--- Final Summary ---")
print(f"The analysis shows that {total_resistant_percentage:.2f}% of offspring should be resistant, only E3ub-wt is an active ligase, only E3ub-wt interacts with Par22, and the insertion adds about {mass_increase_kda:.1f} kDa to the protein mass.")
print("This set of conclusions matches option J.")