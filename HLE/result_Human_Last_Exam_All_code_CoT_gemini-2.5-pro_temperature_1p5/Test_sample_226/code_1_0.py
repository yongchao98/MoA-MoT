# Specimen counts for each family
dixidae_count = 100
simuliidae_count = 1101
blephariceridae_count = 100
rhagionidae_count = 101
tabanidae_count = 201

# Number of prolegs/parapodia per larva for each family
prolegs_dixidae = 4
prolegs_simuliidae = 2
prolegs_blephariceridae = 6
prolegs_rhagionidae = 7
prolegs_tabanidae = 7

# Calculate the subtotal of prolegs for each family
total_dixidae = dixidae_count * prolegs_dixidae
total_simuliidae = simuliidae_count * prolegs_simuliidae
total_blephariceridae = blephariceridae_count * prolegs_blephariceridae
total_rhagionidae = rhagionidae_count * prolegs_rhagionidae
total_tabanidae = tabanidae_count * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blephariceridae + total_rhagionidae + total_tabanidae

# Print the final equation and the answer
print(f"Total Prolegs and Parapodia Calculation:")
print(f"({dixidae_count} Dixidae * {prolegs_dixidae}) + ({simuliidae_count} Simuliidae * {prolegs_simuliidae}) + ({blephariceridae_count} Blephariceridae * {prolegs_blephariceridae}) + ({rhagionidae_count} Rhagionidae * {prolegs_rhagionidae}) + ({tabanidae_count} Tabanidae * {prolegs_tabanidae})")
print(f"= {total_dixidae} + {total_simuliidae} + {total_blephariceridae} + {total_rhagionidae} + {total_tabanidae}")
print(f"= {grand_total}")