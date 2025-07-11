# Step 1: Define the number of specimens for each family.
dixidae_specimens = 100
simuliidae_specimens = 1101
blepharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Step 2: Define the number of prolegs/parapodia per larva for each family.
dixidae_prolegs_per_larva = 6
simuliidae_prolegs_per_larva = 2
blepharoceridae_prolegs_per_larva = 6
rhagionidae_prolegs_per_larva = 0
tabanidae_prolegs_per_larva = 14

# Step 3: Calculate the total number of prolegs for each family.
total_dixidae = dixidae_specimens * dixidae_prolegs_per_larva
total_simuliidae = simuliidae_specimens * simuliidae_prolegs_per_larva
total_blepharoceridae = blepharoceridae_specimens * blepharoceridae_prolegs_per_larva
total_rhagionidae = rhagionidae_specimens * rhagionidae_prolegs_per_larva
total_tabanidae = tabanidae_specimens * tabanidae_prolegs_per_larva

# Step 4: Calculate the grand total.
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Step 5: Print the full calculation and the final single numeric value.
print("Calculating total prolegs and parapodia based on specimen counts:\n")
print(f"Dixidae: {dixidae_specimens} specimens * {dixidae_prolegs_per_larva} prolegs/specimen = {total_dixidae}")
print(f"Simuliidae: {simuliidae_specimens} specimens * {simuliidae_prolegs_per_larva} prolegs/specimen = {total_simuliidae}")
print(f"Blepharoceridae: {blepharoceridae_specimens} specimens * {blepharoceridae_prolegs_per_larva} prolegs/specimen = {total_blepharoceridae}")
print(f"Rhagionidae: {rhagionidae_specimens} specimens * {rhagionidae_prolegs_per_larva} prolegs/specimen = {total_rhagionidae}")
print(f"Tabanidae: {tabanidae_specimens} specimens * {tabanidae_prolegs_per_larva} prolegs/specimen = {total_tabanidae}")
print(f"\nFinal Equation: {total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_rhagionidae} + {total_tabanidae} = {grand_total}")
print(f"\nThe total number of prolegs and parapodia in the dish is {grand_total}.")
