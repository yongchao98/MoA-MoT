# Number of specimens for each family
dixidae_count = 100
simuliidae_count = 1101
belpharoceridae_count = 100
rhagionidae_count = 101
tabanidae_count = 201

# Number of prolegs/parapodia per larva for each family
prolegs_per_dixidae = 4
prolegs_per_simuliidae = 2
prolegs_per_belpharoceridae = 6
prolegs_per_rhagionidae = 0
prolegs_per_tabanidae = 56 # (8 prolegs per segment * 7 segments)

# Calculate the subtotal for each family
dixidae_total = dixidae_count * prolegs_per_dixidae
simuliidae_total = simuliidae_count * prolegs_per_simuliidae
belpharoceridae_total = belpharoceridae_count * prolegs_per_belpharoceridae
rhagionidae_total = rhagionidae_count * prolegs_per_rhagionidae
tabanidae_total = tabanidae_count * prolegs_per_tabanidae

# Calculate the grand total
total_prolegs = dixidae_total + simuliidae_total + belpharoceridae_total + rhagionidae_total + tabanidae_total

# Print the final equation
print("The total number of prolegs and parapodia is calculated as follows:")
print(f"({dixidae_count} * {prolegs_per_dixidae}) + ({simuliidae_count} * {prolegs_per_simuliidae}) + ({belpharoceridae_count} * {prolegs_per_belpharoceridae}) + ({rhagionidae_count} * {prolegs_per_rhagionidae}) + ({tabanidae_count} * {prolegs_per_tabanidae}) = {total_prolegs}")
print("\nTotal number of prolegs and parapodia in the dish:")
print(total_prolegs)
