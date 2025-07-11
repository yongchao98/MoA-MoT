# Constants and given data
prob_self_pollination = 0.05
prob_cross_pollination = 0.95

# The resistant parent is heterozygous (Wt/Ins). We assume resistance is a dominant trait,
# so genotypes Wt/Ins and Ins/Ins are resistant. The general population is Wt/Wt.

# For self-pollination (Wt/Ins x Wt/Ins):
# Offspring: 1/4 Wt/Wt (susceptible), 2/4 Wt/Ins (resistant), 1/4 Ins/Ins (resistant)
resistance_from_selfing = 0.75  # 75%

# For cross-pollination (Wt/Ins x Wt/Wt):
# Offspring: 1/2 Wt/Wt (susceptible), 1/2 Wt/Ins (resistant)
resistance_from_crossing = 0.50  # 50%

# Calculate the total percentage of resistant offspring
total_resistance = (prob_self_pollination * resistance_from_selfing) + (prob_cross_pollination * resistance_from_crossing)
total_resistance_percent = total_resistance * 100

# Calculate the mass increase from the nucleotide insertion
nucleotide_insertion_length = 105
amino_acids_per_codon = 3
average_amino_acid_mass_da = 110 # in Daltons

amino_acid_increase = nucleotide_insertion_length / amino_acids_per_codon
mass_increase_kda = (amino_acid_increase * average_amino_acid_mass_da) / 1000

print("1. Analysis of E3 Ligase Activity and Interaction:")
print("Co-expression data (200 units vs 700) shows E3ub-wt is an active ubiquitin ligase.")
print("Co-expression data (3000 units vs 700) shows E3ub-insert105 is not an active ubiquitin ligase.")
print("Mass spectrometry data (110kDa complex) shows Par22 interacts with E3ub-wt.")
print("Mass spectrometry data (free 50kDa Par22) shows Par22 does NOT interact with E3ub-insert105.\n")

print("2. Calculation of Mass Increase:")
print(f"The 105 nucleotide insertion adds {int(amino_acid_increase)} amino acids.")
print(f"The resulting mass increase is approximately {mass_increase_kda:.1f} kDa (or ~4.0 kDa).\n")

print("3. Calculation of Offspring Resistance:")
print(f"The final equation for the percentage of resistant offspring is:")
print(f"({prob_self_pollination} * {resistance_from_selfing*100}%) + ({prob_cross_pollination} * {resistance_from_crossing*100}%) = {total_resistance_percent:.2f}%")
<<<J>>>