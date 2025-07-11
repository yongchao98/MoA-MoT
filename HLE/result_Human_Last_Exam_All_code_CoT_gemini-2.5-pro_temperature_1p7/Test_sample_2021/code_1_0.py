# Part 1: Calculation of the percentage of resistant offspring

# Given rates of pollination
self_pollination_rate = 0.05  # 5%
cross_pollination_rate = 0.95 # 95%

# Probability of resistant offspring from self-pollination (wt/insert x wt/insert -> 3/4 resistant)
p_resistant_given_self = 0.75

# Probability of resistant offspring from cross-pollination (wt/insert x wt/wt -> 1/2 resistant)
p_resistant_given_cross = 0.50

# Calculate total probability of resistant offspring
total_probability_resistant = (p_resistant_given_self * self_pollination_rate) + (p_resistant_given_cross * cross_pollination_rate)
percentage_resistant = total_probability_resistant * 100

print("--- Offspring Resistance Calculation ---")
print(f"The equation for the total probability of a resistant offspring is:")
print(f"P(Resistant) = (P(Resistant|Self) * P(Self)) + (P(Resistant|Cross) * P(Cross))")
print(f"P(Resistant) = ({p_resistant_given_self} * {self_pollination_rate}) + ({p_resistant_given_cross} * {cross_pollination_rate})")
print(f"P(Resistant) = {p_resistant_given_self * self_pollination_rate} + {p_resistant_given_cross * cross_pollination_rate}")
print(f"P(Resistant) = {total_probability_resistant}")
print(f"Theoretically, {percentage_resistant:.2f}% of the offspring should be drought-resistant.")
print("-" * 35)

# Part 2: Calculation of the mass increase of protein E3ub

# Length of the nucleotide insertion
insertion_length_nt = 105

# Number of amino acids per nucleotide codon
codons_size = 3

# Average mass of an amino acid in Daltons (Da)
avg_amino_acid_mass_da = 110

# Calculate the number of added amino acids
num_amino_acids = insertion_length_nt / codons_size

# Calculate the total mass increase in Daltons
mass_increase_da = num_amino_acids * avg_amino_acid_mass_da

# Convert Daltons to kiloDaltons (kDa)
mass_increase_kda = mass_increase_da / 1000

print("\n--- Protein Mass Increase Calculation ---")
print(f"The insertion has {insertion_length_nt} nucleotides.")
print(f"This corresponds to {insertion_length_nt} / {codons_size} = {int(num_amino_acids)} amino acids.")
print(f"The mass increase is {int(num_amino_acids)} amino acids * {avg_amino_acid_mass_da} Da/amino acid = {int(mass_increase_da)} Da.")
print(f"This is equal to {mass_increase_kda:.2f} kDa, which is approximately 4.0 kDa.")
print("-" * 35)

print("\n--- Final Conclusion ---")
print("Based on the analysis and calculations:")
print("- Theoretically, 51.25% of the offspring should be drought-resistant.")
print("- Only E3ub-wt is an active ubiquitin ligase.")
print("- Par22 cannot interact with E3ub-insert105 (only with E3ub-wt).")
print("- The insertion increases the mass of E3ub by approximately 4.0 kDa.")
print("\nThis corresponds to answer choice J.")