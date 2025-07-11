# Number of specimens collected for each family
dixidae_count = 100
simuliidae_count = 1101
blepharoceridae_count = 100
rhagionidae_count = 101
tabanidae_count = 201

# Number of prolegs/parapodia per larva for each family
# Dixidae larvae have 5 pairs of ventral prolegs.
prolegs_per_dixidae = 10
# Simuliidae larvae have one anterior and one posterior proleg.
prolegs_per_simuliidae = 2
# Blepharoceridae larvae have 6 ventral suctorial discs (modified prolegs).
prolegs_per_blepharoceridae = 6
# Rhagionidae (Vermileoninae) larvae are apodous (no prolegs).
prolegs_per_rhagionidae = 0
# Tabanidae (Tabanus) larvae have rings of 8 prolegs on 7 segments (8*7).
prolegs_per_tabanidae = 56

# Calculate the total number of prolegs for each family
total_dixidae = dixidae_count * prolegs_per_dixidae
total_simuliidae = simuliidae_count * prolegs_per_simuliidae
total_blepharoceridae = blepharoceridae_count * prolegs_per_blepharoceridae
total_rhagionidae = rhagionidae_count * prolegs_per_rhagionidae
total_tabanidae = tabanidae_count * prolegs_per_tabanidae

# Calculate the final total
total_prolegs = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation and the result
print("Calculation:")
print(f"({dixidae_count} * {prolegs_per_dixidae}) + ({simuliidae_count} * {prolegs_per_simuliidae}) + ({blepharoceridae_count} * {prolegs_per_blepharoceridae}) + ({rhagionidae_count} * {prolegs_per_rhagionidae}) + ({tabanidae_count} * {prolegs_per_tabanidae})")
print(f"= {total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_rhagionidae} + {total_tabanidae}")
print(f"= {total_prolegs}")

print("\nFinal Answer:")
print(total_prolegs)