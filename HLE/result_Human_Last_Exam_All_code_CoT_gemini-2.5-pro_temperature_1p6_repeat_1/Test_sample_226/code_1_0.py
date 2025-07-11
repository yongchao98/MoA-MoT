# Number of specimens for each family
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Number of prolegs/parapodia per larva for each family
prolegs_dixidae = 4
prolegs_simuliidae = 2
prolegs_blepharoceridae = 6
prolegs_rhagionidae = 0
prolegs_tabanidae = 14

# Calculate the total number of prolegs for each family
total_dixidae = num_dixidae * prolegs_dixidae
total_simuliidae = num_simuliidae * prolegs_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_blepharoceridae
total_rhagionidae = num_rhagionidae * prolegs_rhagionidae
total_tabanidae = num_tabanidae * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation with the numbers
print(f"({num_dixidae} * {prolegs_dixidae}) + ({num_simuliidae} * {prolegs_simuliidae}) + ({num_blepharoceridae} * {prolegs_blepharoceridae}) + ({num_rhagionidae} * {prolegs_rhagionidae}) + ({num_tabanidae} * {prolegs_tabanidae}) = {grand_total}")

# Print the final answer as a single numeric value
print(f"Total prolegs and parapodia: {grand_total}")