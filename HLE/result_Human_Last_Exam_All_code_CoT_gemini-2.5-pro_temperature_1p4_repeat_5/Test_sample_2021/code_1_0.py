import math

# Step 1: Calculate the mass increase of the E3ub protein.
nucleotide_insertion_length = 105
nucleotides_per_amino_acid = 3
# The average molecular weight of an amino acid is approximately 110 Daltons (Da).
average_amino_acid_mass_da = 110

# Calculate the number of added amino acids.
added_amino_acids = nucleotide_insertion_length / nucleotides_per_amino_acid

# Calculate the total mass increase in kiloDaltons (kDa).
mass_increase_kda = (added_amino_acids * average_amino_acid_mass_da) / 1000

print("--- Analysis of Protein and Genetic Data ---")
print("\nPart 1: Protein Mass Calculation")
print(f"The insertion has {nucleotide_insertion_length} nucleotides.")
print(f"This translates to {int(added_amino_acids)} additional amino acids ({nucleotide_insertion_length} / {nucleotides_per_amino_acid}).")
print(f"The calculated mass increase is approximately {mass_increase_kda:.1f} kDa.")
print("\nConclusion from experiment: The insertion increases the mass of E3ub by approximately 4.0 kDa.")


# Step 2: Calculate the theoretical percentage of drought-resistant offspring.
print("\nPart 2: Offspring Resistance Probability Calculation")

# Pollination rates
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# The parent plant is heterozygous (Ee) and resistant.
# Resistance (e) is dominant, so Ee and ee genotypes are resistant.

# Probability of resistant offspring from self-pollination (Ee x Ee -> 1/4 EE, 2/4 Ee, 1/4 ee)
prob_resistant_from_selfing = 3/4

# Probability of resistant offspring from cross-pollination (Ee x EE -> 1/2 Ee, 1/2 EE)
prob_resistant_from_crossing = 1/2

# Calculate the total probability of resistant offspring.
total_prob_resistant = (prob_resistant_from_selfing * self_pollination_rate) + \
                       (prob_resistant_from_crossing * cross_pollination_rate)

# Convert to percentage.
percentage_resistant_offspring = total_prob_resistant * 100

print("The equation for the final percentage of resistant offspring is:")
# Output each number in the final equation as requested.
print(f"({prob_resistant_from_selfing} * {self_pollination_rate} + {prob_resistant_from_crossing} * {cross_pollination_rate}) * 100")
print(f"\nConclusion from calculation: Theoretically, {percentage_resistant_offspring:.2f}% of the offspring should be drought-resistant.")

print("\n--- Summary of All Conclusions ---")
print("1. Ligase Activity: Only E3ub-wt is an active ubiquitin ligase (based on Par22 degradation).")
print("2. Protein Interaction: Par22 cannot interact with E3ub-insert105 (based on mass spectrometry).")
print(f"3. Mass Increase: The insertion adds ~{mass_increase_kda:.1f} kDa to the protein's mass.")
print(f"4. Offspring Resistance: {percentage_resistant_offspring:.2f}% of offspring should be resistant.")
print("\nThese conclusions align perfectly with answer choice J.")
