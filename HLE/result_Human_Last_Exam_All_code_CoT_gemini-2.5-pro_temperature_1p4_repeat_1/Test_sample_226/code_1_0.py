# Step 1: Define the number of specimens for each family.
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Step 2: Define the number of prolegs/parapodia per larva for each family based on morphology.
# Dixidae larvae have 2 pairs of prolegs = 4 total.
prolegs_per_dixidae = 4
# Simuliidae larvae have 1 anterior proleg and 1 posterior circlet = 2 total.
prolegs_per_simuliidae = 2
# Blepharoceridae larvae have suckers, not prolegs.
prolegs_per_blepharoceridae = 0
# Rhagionidae (Vermileoninae) larvae are apodous (no prolegs).
prolegs_per_rhagionidae = 0
# Tabanidae larvae typically have rings of 4 pairs of prolegs on 7 segments (4*2*7=56).
prolegs_per_tabanidae = 56

# Step 3: Calculate the subtotal of prolegs for each family.
total_dixidae = num_dixidae * prolegs_per_dixidae
total_simuliidae = num_simuliidae * prolegs_per_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_per_blepharoceridae
total_rhagionidae = num_rhagionidae * prolegs_per_rhagionidae
total_tabanidae = num_tabanidae * prolegs_per_tabanidae

# Step 4: Sum the subtotals to find the grand total.
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation with all its components, as requested.
print("Total Prolegs/Parapodia Calculation:")
print(f"({num_dixidae} * {prolegs_per_dixidae}) + ({num_simuliidae} * {prolegs_per_simuliidae}) + ({num_blepharoceridae} * {prolegs_per_blepharoceridae}) + ({num_rhagionidae} * {prolegs_per_rhagionidae}) + ({num_tabanidae} * {prolegs_per_tabanidae})")
print(f"= {total_dixidae} + {total_simuliidae} + {total_blepharoceridae} + {total_rhagionidae} + {total_tabanidae}")
print(f"= {grand_total}")