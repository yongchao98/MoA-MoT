# Step 1: Define specimen counts for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Step 2: Define the number of prolegs/parapodia per larva for each family
prolegs_per_dixidae = 4
prolegs_per_simuliidae = 2
prolegs_per_belpharoceridae = 6
prolegs_per_rhagionidae = 0
prolegs_per_tabanidae = 56 # 4 pairs (8 total) on each of 7 segments = 56

# Step 3: Calculate the total number of prolegs for each family
dixidae_total = dixidae_specimens * prolegs_per_dixidae
simuliidae_total = simuliidae_specimens * prolegs_per_simuliidae
belpharoceridae_total = belpharoceridae_specimens * prolegs_per_belpharoceridae
rhagionidae_total = rhagionidae_specimens * prolegs_per_rhagionidae
tabanidae_total = tabanidae_specimens * prolegs_per_tabanidae

# Step 4: Calculate the grand total
grand_total = dixidae_total + simuliidae_total + belpharoceridae_total + rhagionidae_total + tabanidae_total

# Print the detailed equation and the final result
print(f"Total prolegs and parapodia calculation:")
print(f"({dixidae_specimens} * {prolegs_per_dixidae}) + ({simuliidae_specimens} * {prolegs_per_simuliidae}) + ({belpharoceridae_specimens} * {prolegs_per_belpharoceridae}) + ({rhagionidae_specimens} * {prolegs_per_rhagionidae}) + ({tabanidae_specimens} * {prolegs_per_tabanidae})")
print(f"= {dixidae_total} + {simuliidae_total} + {belpharoceridae_total} + {rhagionidae_total} + {tabanidae_total}")
print(f"= {grand_total}")