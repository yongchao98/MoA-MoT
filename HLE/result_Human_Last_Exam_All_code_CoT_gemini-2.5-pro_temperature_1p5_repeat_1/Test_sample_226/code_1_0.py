# Specimen counts from the collection data
num_dixidae = 100
num_simuliidae = 1101
num_belpharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Number of prolegs/parapodia per larva for each family based on morphology
prolegs_dixidae = 4
prolegs_simuliidae = 2
prolegs_belpharoceridae = 6
prolegs_rhagionidae = 7
prolegs_tabanidae = 56 # (4 pairs per segment * 7 segments)

# Calculate the subtotal of prolegs for each family
total_dixidae = num_dixidae * prolegs_dixidae
total_simuliidae = num_simuliidae * prolegs_simuliidae
total_belpharoceridae = num_belpharoceridae * prolegs_belpharoceridae
total_rhagionidae = num_rhagionidae * prolegs_rhagionidae
total_tabanidae = num_tabanidae * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae

# Print the full equation showing each number used in the calculation
print(f"The calculation for the total number of prolegs and parapodia is:")
print(f"({num_dixidae} * {prolegs_dixidae}) + ({num_simuliidae} * {prolegs_simuliidae}) + ({num_belpharoceridae} * {prolegs_belpharoceridae}) + ({num_rhagionidae} * {prolegs_rhagionidae}) + ({num_tabanidae} * {prolegs_tabanidae}) = {grand_total}")

# Print the single numeric value of the answer
print("\nTotal number of prolegs and parapodia:")
print(grand_total)