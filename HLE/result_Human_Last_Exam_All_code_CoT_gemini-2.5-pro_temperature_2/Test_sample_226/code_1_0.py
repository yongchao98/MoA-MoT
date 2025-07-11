# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
blepharoceridae_specimens = 100
vermileoninae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per larva for each family
prolegs_per_dixidae = 6
prolegs_per_simuliidae = 2
prolegs_per_blepharoceridae = 6
prolegs_per_vermileoninae = 0
prolegs_per_tabanidae = 42 # 3 pairs on each of 7 abdominal segments

# Calculate the subtotal for each family
total_dixidae = dixidae_specimens * prolegs_per_dixidae
total_simuliidae = simuliidae_specimens * prolegs_per_simuliidae
total_blepharoceridae = blepharoceridae_specimens * prolegs_per_blepharoceridae
total_vermileoninae = vermileoninae_specimens * prolegs_per_vermileoninae
total_tabanidae = tabanidae_specimens * prolegs_per_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_vermileoninae + total_tabanidae

# Print the breakdown of the calculation
print("Calculation Breakdown:")
print(f"Dixidae: {dixidae_specimens} specimens * {prolegs_per_dixidae} prolegs/specimen = {total_dixidae} prolegs")
print(f"Simuliidae: {simuliidae_specimens} specimens * {prolegs_per_simuliidae} prolegs/specimen = {total_simuliidae} prolegs")
print(f"Blepharoceridae: {blepharoceridae_specimens} specimens * {prolegs_per_blepharoceridae} prolegs/specimen = {total_blepharoceridae} prolegs")
print(f"Vermileoninae: {vermileoninae_specimens} specimens * {prolegs_per_vermileoninae} prolegs/specimen = {total_vermileoninae} prolegs")
print(f"Tabanidae: {tabanidae_specimens} specimens * {prolegs_per_tabanidae} prolegs/specimen = {total_tabanidae} prolegs")
print("\nFinal Equation:")
print(f"{total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_vermileoninae} + {total_tabanidae} = {grand_total}")

print("\nTotal number of prolegs and parapodia:")
print(grand_total)
