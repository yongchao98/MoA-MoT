# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
belpharoceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per larva for each family
dixidae_prolegs_per_larva = 4
simuliidae_prolegs_per_larva = 2
belpharoceridae_prolegs_per_larva = 6
rhagionidae_prolegs_per_larva = 0
tabanidae_prolegs_per_larva = 7

# Calculate the total number of prolegs for each family
dixidae_total = dixidae_specimens * dixidae_prolegs_per_larva
simuliidae_total = simuliidae_specimens * simuliidae_prolegs_per_larva
belpharoceridae_total = belpharoceridae_specimens * belpharoceridae_prolegs_per_larva
rhagionidae_total = rhagionidae_specimens * rhagionidae_prolegs_per_larva
tabanidae_total = tabanidae_specimens * tabanidae_prolegs_per_larva

# Calculate the grand total
total_prolegs = dixidae_total + simuliidae_total + belpharoceridae_total + rhagionidae_total + tabanidae_total

# Print the breakdown of the calculation and the final answer
print("Calculation Breakdown:")
print(f"Dixidae: {dixidae_specimens} specimens * {dixidae_prolegs_per_larva} prolegs/specimen = {dixidae_total} prolegs")
print(f"Simuliidae: {simuliidae_specimens} specimens * {simuliidae_prolegs_per_larva} prolegs/specimen = {simuliidae_total} prolegs")
print(f"Belpharoceridae: {belpharoceridae_specimens} specimens * {belpharoceridae_prolegs_per_larva} prolegs/specimen = {belpharoceridae_total} prolegs")
print(f"Rhagionidae: {rhagionidae_specimens} specimens * {rhagionidae_prolegs_per_larva} prolegs/specimen = {rhagionidae_total} prolegs")
print(f"Tabanidae: {tabanidae_specimens} specimens * {tabanidae_prolegs_per_larva} prolegs/specimen = {tabanidae_total} prolegs")
print("\nFinal Equation:")
print(f"({dixidae_specimens} * {dixidae_prolegs_per_larva}) + ({simuliidae_specimens} * {simuliidae_prolegs_per_larva}) + ({belpharoceridae_specimens} * {belpharoceridae_prolegs_per_larva}) + ({rhagionidae_specimens} * {rhagionidae_prolegs_per_larva}) + ({tabanidae_specimens} * {tabanidae_prolegs_per_larva}) = {total_prolegs}")
print(f"\nTotal number of prolegs and parapodia: {total_prolegs}")