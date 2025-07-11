# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
blepharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs per larva for each family
prolegs_dixidae = 4
prolegs_simuliidae = 2
prolegs_blepharoceridae = 0
prolegs_rhagionidae = 0
prolegs_tabanidae = 56 # (4 pairs or 8 prolegs per segment * 7 segments)

# Calculate the total prolegs for each family
total_dixidae = dixidae_specimens * prolegs_dixidae
total_simuliidae = simuliidae_specimens * prolegs_simuliidae
total_blepharoceridae = blepharoceridae_specimens * prolegs_blepharoceridae
total_rhagionidae = rhagionidae_specimens * prolegs_rhagionidae
total_tabanidae = tabanidae_specimens * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the breakdown of the calculation
print("Calculation Breakdown:")
print(f"Dixidae: {dixidae_specimens} specimens * {prolegs_dixidae} prolegs/specimen = {total_dixidae}")
print(f"Simuliidae: {simuliidae_specimens} specimens * {prolegs_simuliidae} prolegs/specimen = {total_simuliidae}")
print(f"Blepharoceridae: {blepharoceridae_specimens} specimens * {prolegs_blepharoceridae} prolegs/specimen = {total_blepharoceridae}")
print(f"Rhagionidae: {rhagionidae_specimens} specimens * {prolegs_rhagionidae} prolegs/specimen = {total_rhagionidae}")
print(f"Tabanidae: {tabanidae_specimens} specimens * {prolegs_tabanidae} prolegs/specimen = {total_tabanidae}")
print("\nFinal Equation:")
print(f"({dixidae_specimens} * {prolegs_dixidae}) + ({simuliidae_specimens} * {prolegs_simuliidae}) + ({blepharoceridae_specimens} * {prolegs_blepharoceridae}) + ({rhagionidae_specimens} * {prolegs_rhagionidae}) + ({tabanidae_specimens} * {prolegs_tabanidae}) = {grand_total}")
print("\nTotal prolegs and parapodia in the dish:")
print(grand_total)
