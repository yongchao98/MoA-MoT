# Number of specimens for each family
num_dixidae = 100
num_simuliidae = 1101
num_belpharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Number of prolegs/parapodia per larva for each family
prolegs_per_dixidae = 4
prolegs_per_simuliidae = 2
prolegs_per_belpharoceridae = 0 # They have suckers, not prolegs
prolegs_per_rhagionidae = 0 # They are apodous (legless)
prolegs_per_tabanidae = 56 # 4 pairs of prolegs on each of the first 7 abdominal segments (4*2*7)

# Calculate the total number of prolegs for each family
total_dixidae = num_dixidae * prolegs_per_dixidae
total_simuliidae = num_simuliidae * prolegs_per_simuliidae
total_belpharoceridae = num_belpharoceridae * prolegs_per_belpharoceridae
total_rhagionidae = num_rhagionidae * prolegs_per_rhagionidae
total_tabanidae = num_tabanidae * prolegs_per_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae

# Print the equation and the final answer
print("Calculating the total number of prolegs and parapodia:")
print(f"Dixidae: {num_dixidae} specimens * {prolegs_per_dixidae} prolegs/specimen = {total_dixidae}")
print(f"Simuliidae: {num_simuliidae} specimens * {prolegs_per_simuliidae} prolegs/specimen = {total_simuliidae}")
print(f"Belpharoceridae: {num_belpharoceridae} specimens * {prolegs_per_belpharoceridae} prolegs/specimen = {total_belpharoceridae}")
print(f"Rhagionidae: {num_rhagionidae} specimens * {prolegs_per_rhagionidae} prolegs/specimen = {total_rhagionidae}")
print(f"Tabanidae: {num_tabanidae} specimens * {prolegs_per_tabanidae} prolegs/specimen = {total_tabanidae}")
print("-" * 30)
print(f"Total = {total_dixidae} + {total_simuliidae} + {total_belpharoceridae} + {total_rhagionidae} + {total_tabanidae} = {grand_total}")
print(f"\nTotal number of prolegs and parapodia in the dish: {grand_total}")
