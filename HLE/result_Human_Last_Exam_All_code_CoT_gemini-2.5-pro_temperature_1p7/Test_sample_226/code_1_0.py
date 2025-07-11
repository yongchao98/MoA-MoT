# Number of specimens for each family
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Number of prolegs/parapodia per larval specimen
# Dixidae larvae have 2 pairs of prolegs.
prolegs_per_dixidae = 4
# Simuliidae larvae have 1 anterior proleg and 1 posterior anal proleg.
prolegs_per_simuliidae = 2
# Blepharoceridae larvae have 6 ventral suckers.
prolegs_per_blepharoceridae = 6
# Rhagionidae (Vermileoninae) larvae are apodous (legless).
prolegs_per_rhagionidae = 0
# Tabanidae larvae have prolegs (creeping welts) on the first 7 abdominal segments.
# We are counting the ventral pair on each of these segments (7 pairs = 14 prolegs).
prolegs_per_tabanidae = 14

# Calculate the total prolegs for each family
total_dixidae = num_dixidae * prolegs_per_dixidae
total_simuliidae = num_simuliidae * prolegs_per_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_per_blepharoceridae
total_rhagionidae = num_rhagionidae * prolegs_per_rhagionidae
total_tabanidae = num_tabanidae * prolegs_per_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation with each component number
print("Final Equation:")
print(f"{total_dixidae} (from Dixidae) + {total_simuliidae} (from Simuliidae) + {total_blepharoceridae} (from Blepharoceridae) + {total_rhagionidae} (from Rhagionidae) + {total_tabanidae} (from Tabanidae) = {grand_total}")
print("\nTotal number of prolegs and parapodia:")
print(grand_total)