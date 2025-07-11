# Number of specimens for each family
dixidae_specimens = 100
simuliidae_specimens = 1101
blephariceridae_specimens = 100
rhagionidae_specimens = 101
tabanidae_specimens = 201

# Number of prolegs/parapodia per specimen for each family
prolegs_per_dixidae = 6       # 3 pairs of prolegs
prolegs_per_simuliidae = 2    # 1 anterior and 1 posterior proleg
prolegs_per_blephariceridae = 12  # 6 pairs of prolegs
prolegs_per_rhagionidae = 0   # Apodous (no prolegs)
prolegs_per_tabanidae = 42    # 3 pairs on each of 7 segments

# Calculate the subtotal for each family
dixidae_total = dixidae_specimens * prolegs_per_dixidae
simuliidae_total = simuliidae_specimens * prolegs_per_simuliidae
blephariceridae_total = blephariceridae_specimens * prolegs_per_blephariceridae
rhagionidae_total = rhagionidae_specimens * prolegs_per_rhagionidae
tabanidae_total = tabanidae_specimens * prolegs_per_tabanidae

# Calculate the grand total
grand_total = dixidae_total + simuliidae_total + blephariceridae_total + rhagionidae_total + tabanidae_total

# Print the breakdown of the calculation
print(f"Dixidae: {dixidae_specimens} specimens * {prolegs_per_dixidae} prolegs/specimen = {dixidae_total} prolegs")
print(f"Simuliidae: {simuliidae_specimens} specimens * {prolegs_per_simuliidae} prolegs/specimen = {simuliidae_total} prolegs")
print(f"Blephariceridae: {blephariceridae_specimens} specimens * {prolegs_per_blephariceridae} prolegs/specimen = {blephariceridae_total} prolegs")
print(f"Rhagionidae (Vermileoninae): {rhagionidae_specimens} specimens * {prolegs_per_rhagionidae} prolegs/specimen = {rhagionidae_total} prolegs")
print(f"Tabanidae (Tabanus): {tabanidae_specimens} specimens * {prolegs_per_tabanidae} prolegs/specimen = {tabanidae_total} prolegs")
print("\nCalculating the total sum:")
print(f"{dixidae_total} + {simuliidae_total} + {blephariceridae_total} + {rhagionidae_total} + {tabanidae_total} = {grand_total}")
print(f"\nTotal number of prolegs and parapodia in the dish: {grand_total}")
