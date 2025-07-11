# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per specimen for each family
prolegs_per_dixidae = 4
prolegs_per_simuliidae = 2
prolegs_per_belpharoceridae = 0
prolegs_per_rhagionidae = 7
prolegs_per_tabanidae = 56

# Calculate the total for each family
dixidae_total = dixidae_specimens * prolegs_per_dixidae
simuliidae_total = simuliidae_specimens * prolegs_per_simuliidae
belpharoceridae_total = belpharoceridae_specimens * prolegs_per_belpharoceridae
rhagionidae_total = rhagionidae_specimens * prolegs_per_rhagionidae
tabanidae_total = tabanidae_specimens * prolegs_per_tabanidae

# Calculate the grand total
grand_total = dixidae_total + simuliidae_total + belpharoceridae_total + rhagionidae_total + tabanidae_total

# Print the final equation as requested
print(f"({dixidae_specimens} * {prolegs_per_dixidae}) + ({simuliidae_specimens} * {prolegs_per_simuliidae}) + ({belpharoceridae_specimens} * {prolegs_per_belpharoceridae}) + ({rhagionidae_specimens} * {prolegs_per_rhagionidae}) + ({tabanidae_specimens} * {prolegs_per_tabanidae}) = {grand_total}")

# Print the final numeric value
print(f"Total number of prolegs and parapodia: {grand_total}")
<<<14565>>>