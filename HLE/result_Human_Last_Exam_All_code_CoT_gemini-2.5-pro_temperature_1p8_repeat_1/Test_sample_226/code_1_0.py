# Number of specimens for each family
dixidae_count = 100
simuliidae_count = 1101
belpharoceridae_count = 100
rhagionidae_count = 101
tabanidae_count = 201

# Number of prolegs/parapodia per larva for each family
prolegs_dixidae = 6
prolegs_simuliidae = 2
prolegs_belpharoceridae = 6
prolegs_rhagionidae = 0
prolegs_tabanidae = 56

# Calculate the subtotal for each family
total_dixidae = dixidae_count * prolegs_dixidae
total_simuliidae = simuliidae_count * prolegs_simuliidae
total_belpharoceridae = belpharoceridae_count * prolegs_belpharoceridae
total_rhagionidae = rhagionidae_count * prolegs_rhagionidae
total_tabanidae = tabanidae_count * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation and the answer
print(f"({dixidae_count} * {prolegs_dixidae}) + ({simuliidae_count} * {prolegs_simuliidae}) + ({belpharoceridae_count} * {prolegs_belpharoceridae}) + ({rhagionidae_count} * {prolegs_rhagionidae}) + ({tabanidae_count} * {prolegs_tabanidae}) = {grand_total}")
print("\nTotal number of prolegs and parapodia:")
print(grand_total)
