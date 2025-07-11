# Number of specimens for each family
dixidae_spp_count = 100
simuliidae_spp_count = 1101
belpharoceridae_spp_count = 100
rhagionidae_spp_count = 101
tabanidae_spp_count = 201

# Number of prolegs/parapodia per larva for each family
prolegs_dixidae = 4  # Two pairs of ventral prolegs
prolegs_simuliidae = 2  # One anterior proleg and one posterior anal proleg
prolegs_belpharoceridae = 0  # Have suckers, not prolegs
prolegs_rhagionidae = 0  # Apodous (legless)
prolegs_tabanidae = 7  # Creeping welts on the first 7 abdominal segments

# Calculate the total for each family
total_dixidae = dixidae_spp_count * prolegs_dixidae
total_simuliidae = simuliidae_spp_count * prolegs_simuliidae
total_belpharoceridae = belpharoceridae_spp_count * prolegs_belpharoceridae
total_rhagionidae = rhagionidae_spp_count * prolegs_rhagionidae
total_tabanidae = tabanidae_spp_count * prolegs_tabanidae

# Calculate the grand total
grand_total = total_dixidae + total_simuliidae + total_belpharoceridae + total_rhagionidae + total_tabanidae

# Print the full equation and the final result
print(f"({dixidae_spp_count} * {prolegs_dixidae}) + "
      f"({simuliidae_spp_count} * {prolegs_simuliidae}) + "
      f"({belpharoceridae_spp_count} * {prolegs_belpharoceridae}) + "
      f"({rhagionidae_spp_count} * {prolegs_rhagionidae}) + "
      f"({tabanidae_spp_count} * {prolegs_tabanidae}) = {grand_total}")
print(f"Total prolegs and parapodia = {grand_total}")
