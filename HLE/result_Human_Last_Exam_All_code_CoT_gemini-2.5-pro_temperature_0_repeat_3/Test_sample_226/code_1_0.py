# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
blepharoceridae_specimens = 100
vermileoninae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per specimen for each family
dixidae_prolegs = 4
simuliidae_prolegs = 2
blepharoceridae_prolegs = 6
vermileoninae_prolegs = 0
tabanidae_prolegs = 7

# Calculate the total for each family
total_dixidae = dixidae_specimens * dixidae_prolegs
total_simuliidae = simuliidae_specimens * simuliidae_prolegs
total_blepharoceridae = blepharoceridae_specimens * blepharoceridae_prolegs
total_vermileoninae = vermileoninae_specimens * vermileoninae_prolegs
total_tabanidae = tabanidae_specimens * tabanidae_prolegs

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_vermileoninae + total_tabanidae

# Print the final equation and the result
print(f"Total prolegs and parapodia = ({dixidae_specimens} * {dixidae_prolegs}) + ({simuliidae_specimens} * {simuliidae_prolegs}) + ({blepharoceridae_specimens} * {blepharoceridae_prolegs}) + ({vermileoninae_specimens} * {vermileoninae_prolegs}) + ({tabanidae_specimens} * {tabanidae_prolegs})")
print(f"Total = {total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_vermileoninae} + {total_tabanidae}")
print(f"The total number of prolegs and parapodia in the dish is: {grand_total}")
