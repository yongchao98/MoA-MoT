# Number of specimens for each family
num_dixidae = 100
num_simuliidae = 1101
num_blepharoceridae = 100
num_rhagionidae = 101
num_tabanidae = 201

# Number of prolegs and parapodia per larva for each family based on typical morphology
prolegs_dixidae = 10        # 2 pairs of prolegs (4) + 3 pairs of parapodia (6)
prolegs_simuliidae = 2      # 1 anterior proleg + 1 posterior circlet
prolegs_blepharoceridae = 6   # 6 median ventral suckers
prolegs_rhagionidae = 7     # 7 ventral creeping welts
prolegs_tabanidae = 56      # 4 pairs of prolegs on each of the first 7 abdominal segments (4*2*7)

# Calculate the total for each family
total_dixidae = num_dixidae * prolegs_dixidae
total_simuliidae = num_simuliidae * prolegs_simuliidae
total_blepharoceridae = num_blepharoceridae * prolegs_blepharoceridae
total_rhagionidae = num_rhagionidae * prolegs_rhagionidae
total_tabanidae = num_tabanidae * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_blepharoceridae + total_rhagionidae + total_tabanidae

# Print the final equation and the answer
print("To find the total number of prolegs and parapodia, we sum the totals for each family.")
print("The calculation is as follows:")
print(f"({num_dixidae} Dixidae * {prolegs_dixidae}) + ({num_simuliidae} Simuliidae * {prolegs_simuliidae}) + ({num_blepharoceridae} Blepharoceridae * {prolegs_blepharoceridae}) + ({num_rhagionidae} Rhagionidae * {prolegs_rhagionidae}) + ({num_tabanidae} Tabanidae * {prolegs_tabanidae})")
print(f"= ({total_dixidae}) + ({total_simuliidae}) + ({total_blepharoceridae}) + ({total_rhagionidae}) + ({total_tabanidae})")
print(f"= {grand_total}")
print(f"\nThe total number of prolegs and parapodia in the dish is: {grand_total}")