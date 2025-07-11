# Step 1: Define the number of specimens for each family.
num_dixidae = 100
num_simuliidae = 1101
num_blephariceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Step 2: Define the number of prolegs/parapodia per larva for each family.
# Dixidae larvae have 2 pairs of prolegs.
prolegs_per_dixidae = 4
# Simuliidae larvae have an anterior proleg and a posterior circlet.
prolegs_per_simuliidae = 2
# Blephariceridae larvae have 6 ventral suckers (modified prolegs).
prolegs_per_blephariceridae = 6
# Rhagionidae (Vermileoninae) larvae are apodous (have no prolegs).
prolegs_per_rhagionidae = 0
# Tabanidae larvae have rings of prolegs (8 per ring) on 7 segments.
prolegs_per_tabanidae = 56

# Step 3: Calculate the subtotal of prolegs for each family.
total_dixidae = num_dixidae * prolegs_per_dixidae
total_simuliidae = num_simuliidae * prolegs_per_simuliidae
total_blephariceridae = num_blephariceridae * prolegs_per_blephariceridae
total_rhagionidae = num_rhagionidae * prolegs_per_rhagionidae
total_tabanidae = num_tabanidae * prolegs_per_tabanidae

# Step 4: Calculate the grand total.
grand_total = total_dixidae + total_simuliidae + total_blephariceridae + total_rhagionidae + total_tabanidae

# Step 5: Print the full equation and the final result.
print("Total prolegs and parapodia calculation:")
print(f"({num_dixidae} Dixidae * {prolegs_per_dixidae} prolegs) + ({num_simuliidae} Simuliidae * {prolegs_per_simuliidae} prolegs) + ({num_blephariceridae} Blephariceridae * {prolegs_per_blephariceridae} prolegs) + ({num_rhagionidae} Rhagionidae * {prolegs_per_rhagionidae} prolegs) + ({num_tabanidae} Tabanidae * {prolegs_per_tabanidae} prolegs)")
print(f"= {total_dixidae} + {total_simuliidae} + {total_blephariceridae} + {total_rhagionidae} + {total_tabanidae}")
print(f"= {grand_total}")

print("\nTotal number of prolegs and parapodia in the dish:")
print(grand_total)