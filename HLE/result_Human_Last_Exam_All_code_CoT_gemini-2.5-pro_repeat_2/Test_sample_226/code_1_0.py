# Step 1: Define the number of specimens for each family.
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Step 2: Define the number of prolegs/parapodia for a single larva of each family.
prolegs_per_dixidae = 4       # 2 pairs of prolegs
prolegs_per_simuliidae = 1    # 1 single prothoracic proleg
prolegs_per_belpharoceridae = 0 # Have suckers, not prolegs
prolegs_per_rhagionidae = 0     # Apodous (legless)
prolegs_per_tabanidae = 42      # 6 prolegs per segment on 7 segments

# Step 3: Calculate the subtotal for each family.
dixidae_total = dixidae_specimens * prolegs_per_dixidae
simuliidae_total = simuliidae_specimens * prolegs_per_simuliidae
belpharoceridae_total = belpharoceridae_specimens * prolegs_per_belpharoceridae
rhagionidae_total = rhagionidae_specimens * prolegs_per_rhagionidae
tabanidae_total = tabanidae_specimens * prolegs_per_tabanidae

# Step 4: Calculate the grand total.
grand_total = dixidae_total + simuliidae_total + belpharoceridae_total + rhagionidae_total + tabanidae_total

# Print the final equation with all numbers.
print("Final Equation:")
print(f"({dixidae_specimens} * {prolegs_per_dixidae}) + ({simuliidae_specimens} * {prolegs_per_simuliidae}) + ({belpharoceridae_specimens} * {prolegs_per_belpharoceridae}) + ({rhagionidae_specimens} * {prolegs_per_rhagionidae}) + ({tabanidae_specimens} * {prolegs_per_tabanidae}) = {grand_total}")

print("\nTotal number of prolegs and parapodia:")
print(grand_total)