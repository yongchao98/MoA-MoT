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
prolegs_per_tabanidae = 56

# Calculate the total for each family
total_dixidae = dixidae_specimens * prolegs_per_dixidae
total_simuliidae = simuliidae_specimens * prolegs_per_simuliidae
total_blepharoceridae = blepharoceridae_specimens * prolegs_per_blepharoceridae
total_vermileoninae = vermileoninae_specimens * prolegs_per_vermileoninae
total_tabanidae = tabanidae_specimens * prolegs_per_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_vermileoninae + total_tabanidae

# Print the full equation and the final answer
print("Calculation of total prolegs and parapodia:")
print(f"({dixidae_specimens} Dixidae * {prolegs_per_dixidae} prolegs) + "
      f"({simuliidae_specimens} Simuliidae * {prolegs_per_simuliidae} prolegs) + "
      f"({blepharoceridae_specimens} Belpharoceridae * {prolegs_per_blepharoceridae} prolegs) + "
      f"({vermileoninae_specimens} Rhagionidae * {prolegs_per_vermileoninae} prolegs) + "
      f"({tabanidae_specimens} Tabanidae * {prolegs_per_tabanidae} prolegs)")

print(f"= {total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_vermileoninae} + {total_tabanidae}")
print(f"= {grand_total}")

print("\nTotal number of prolegs and parapodia in the dish:")
print(grand_total)