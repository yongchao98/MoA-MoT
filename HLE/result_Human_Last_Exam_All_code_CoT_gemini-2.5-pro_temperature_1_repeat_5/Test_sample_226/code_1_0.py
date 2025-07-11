# --- Number of specimens for each family ---
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# --- Number of prolegs/parapodia per larva for each family ---
# Dixidae larvae have 2 pairs of ventral prolegs.
prolegs_per_dixidae = 4
# Simuliidae larvae have 1 anterior proleg and 1 posterior circlet.
prolegs_per_simuliidae = 2
# Blephariceridae larvae have suctorial discs, not prolegs.
prolegs_per_blepharoceridae = 0
# Vermileoninae larvae are vermiform (worm-like) and lack prolegs.
prolegs_per_rhagionidae = 0
# Tabanidae larvae have rings of prolegs on 7 abdominal segments (8 prolegs per ring).
prolegs_per_tabanidae = 56

# --- Calculate subtotals for each family ---
total_dixidae = num_dixidae * prolegs_per_dixidae
total_simuliidae = num_simuliidae * prolegs_per_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_per_blepharoceridae
total_rhagionidae = num_rhagionidae * prolegs_per_rhagionidae
total_tabanidae = num_tabanidae * prolegs_per_tabanidae

# --- Calculate the grand total ---
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# --- Print the equation and the final answer ---
print("Calculation of total prolegs and parapodia:")
print(f"({num_dixidae} Dixidae * {prolegs_per_dixidae} prolegs) + ({num_simuliidae} Simuliidae * {prolegs_per_simuliidae} prolegs) + ({num_blepharoceridae} Blepharoceridae * {prolegs_per_blepharoceridae} prolegs) + ({num_rhagionidae} Rhagionidae * {prolegs_per_rhagionidae} prolegs) + ({num_tabanidae} Tabanidae * {prolegs_per_tabanidae} prolegs) = {grand_total}")
print("\nTotal number of prolegs and parapodia in the dish:")
print(grand_total)
