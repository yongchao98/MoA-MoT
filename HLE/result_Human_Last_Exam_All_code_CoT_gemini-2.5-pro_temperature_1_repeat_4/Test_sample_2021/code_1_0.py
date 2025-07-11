# Define pollination rates
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# Define the probability of getting a resistant offspring from each type of cross.
# This assumes that both heterozygous (wt/insert) and homozygous (insert/insert) are resistant.
# For a self-cross (wt/insert x wt/insert), offspring are 1/4 wt/wt, 2/4 wt/insert, 1/4 insert/insert.
# So, the resistant fraction is 2/4 + 1/4 = 3/4 = 0.75.
prob_resistant_from_selfing = 0.75

# For a cross with a wild-type plant (wt/insert x wt/wt), offspring are 1/2 wt/wt, 1/2 wt/insert.
# So, the resistant fraction is 1/2 = 0.5.
prob_resistant_from_crossing = 0.5

# Calculate the total probability of an offspring being resistant
total_prob_resistant = (self_pollination_rate * prob_resistant_from_selfing) + (cross_pollination_rate * prob_resistant_from_crossing)

# Convert to percentage
percentage_resistant = total_prob_resistant * 100

print("Calculation for the percentage of drought-resistant offspring:")
print(f"({self_pollination_rate} * {prob_resistant_from_selfing}) + ({cross_pollination_rate} * {prob_resistant_from_crossing}) = {total_prob_resistant}")
print(f"Theoretically, {percentage_resistant:.2f}% of the offspring should be drought-resistant.\n")

# Calculation for the protein mass increase
insertion_nucleotides = 105
nucleotides_per_amino_acid = 3
avg_mass_per_amino_acid_da = 110 # in Daltons

# Calculate number of amino acids added
added_amino_acids = insertion_nucleotides / nucleotides_per_amino_acid

# Calculate mass increase in Daltons and kiloDaltons
mass_increase_da = added_amino_acids * avg_mass_per_amino_acid_da
mass_increase_kda = mass_increase_da / 1000

print("Calculation for the approximate mass increase of the protein:")
print(f"Number of added amino acids = {insertion_nucleotides} nucleotides / {nucleotides_per_amino_acid} = {int(added_amino_acids)}")
print(f"Mass increase ≈ {int(added_amino_acids)} amino acids * {avg_mass_per_amino_acid_da} Da = {int(mass_increase_da)} Da")
print(f"This is approximately {mass_increase_kda:.2f} kDa, which rounds to 4.0 kDa.")
