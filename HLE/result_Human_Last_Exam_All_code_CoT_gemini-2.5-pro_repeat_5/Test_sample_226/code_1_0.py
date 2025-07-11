# Number of specimens for each family
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_vermileoninae = 101
num_tabanidae = 201

# Number of prolegs/parapodia per larva for each family
prolegs_dixidae = 4
prolegs_simuliidae = 2
prolegs_blepharoceridae = 6
prolegs_vermileoninae = 0
prolegs_tabanidae = 56 # 8 prolegs per segment * 7 segments

# Calculate the subtotal for each family
total_dixidae = num_dixidae * prolegs_dixidae
total_simuliidae = num_simuliidae * prolegs_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_blepharoceridae
total_vermileoninae = num_vermileoninae * prolegs_vermileoninae
total_tabanidae = num_tabanidae * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_vermileoninae + total_tabanidae

# Print the equation
print("Calculation Breakdown:")
print(f"Dixidae: {num_dixidae} specimens * {prolegs_dixidae} prolegs/specimen = {total_dixidae} prolegs")
print(f"Simuliidae: {num_simuliidae} specimens * {prolegs_simuliidae} prolegs/specimen = {total_simuliidae} prolegs")
print(f"Blepharoceridae: {num_blepharoceridae} specimens * {prolegs_blepharoceridae} prolegs/specimen = {total_blepharoceridae} prolegs")
print(f"Rhagionidae (Vermileoninae): {num_vermileoninae} specimens * {prolegs_vermileoninae} prolegs/specimen = {total_vermileoninae} prolegs")
print(f"Tabanidae: {num_tabanidae} specimens * {prolegs_tabanidae} prolegs/specimen = {total_tabanidae} prolegs")
print("\nFinal Equation:")
print(f"{total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_vermileoninae} + {total_tabanidae} = {grand_total}")

print("\nTotal number of prolegs and parapodia:")
print(grand_total)