# --- Problem Variables ---
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
insertion_length_nt = 105
nt_per_amino_acid = 3
avg_mass_per_amino_acid_da = 110

# --- Part 1: Calculate Mass Increase ---
# Calculate the number of amino acids added by the insertion
added_amino_acids = insertion_length_nt / nt_per_amino_acid
# Calculate the total mass increase in Daltons
mass_increase_da = added_amino_acids * avg_mass_per_amino_acid_da
# Convert the mass to kiloDaltons (kDa)
mass_increase_kda = mass_increase_da / 1000

# --- Part 2: Calculate Frequency of Resistant Offspring ---
# From the problem, we deduce that heterozygous (Wt/Ins) and homozygous mutant (Ins/Ins)
# genotypes are resistant, while homozygous wild-type (Wt/Wt) is not.
# Resistant fraction from self-pollination (Wt/Ins x Wt/Ins -> 1/4 Wt/Wt, 2/4 Wt/Ins, 1/4 Ins/Ins)
# Resistant genotypes are Wt/Ins and Ins/Ins, so the fraction is 2/4 + 1/4 = 0.75
resistant_fraction_from_selfing = 0.75
# Resistant fraction from cross-pollination with wild-type (Wt/Ins x Wt/Wt -> 1/2 Wt/Wt, 1/2 Wt/Ins)
# The only resistant genotype is Wt/Ins, so the fraction is 1/2 = 0.5
resistant_fraction_from_crossing = 0.5

# Calculate the total theoretical percentage of resistant offspring
total_resistant_percentage = (self_pollination_rate * resistant_fraction_from_selfing +
                              cross_pollination_rate * resistant_fraction_from_crossing) * 100

# --- Part 3: Formulate the final conclusion ---
# The analysis shows that only E3ub-wt is an active ligase and can interact with Par22.
# E3ub-insert105 is inactive and does not interact.
final_conclusion = (
    f"Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant. "
    f"Only E3ub-wt is an active ubiquitin ligase. "
    f"Par22 cannot interact with E3ub-insert105. "
    f"The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f}kDA."
)

print(final_conclusion)