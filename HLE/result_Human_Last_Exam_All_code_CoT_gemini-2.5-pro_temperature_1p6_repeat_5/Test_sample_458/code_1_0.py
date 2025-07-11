# The task is to identify dragonfly species with reduced pterostigmata based on their ecology.
# A pterostigma helps stabilize the wing during high-speed flight and gliding.
# Therefore, species that are weaker fliers or have specialized flight patterns that do not involve high-speed gliding are likely to have reduced pterostigmata.

# Analysis of the species:
# 1. Didymops transversa (Stream Cruiser): An active, cruising flier. Likely has a normal pterostigma.
# 2. Urothemis edwarsi (Indigo Dropwing): A percher, but not noted for especially weak flight.
# 3. Macrodiplax balteata (Marl Pennant): A strong migrant. Expects a well-developed pterostigma.
# 4. Pantala flavescens (Wandering Glider): The quintessential migrant/glider. Expects a well-developed pterostigma.
# 5. Orthetrum cancellatum (Black-tailed Skimmer): A strong, robust percher.
# 6. Libellula quadrimaculata (Four-spotted Chaser): A known migrant. Expects a well-developed pterostigma.
# 7. Libellula pulchela (Twelve-spotted Skimmer): A strong, territorial percher.
# 8. Sympetrum corruptum (Variegated Meadowhawk): A highly migratory species. Expects a well-developed pterostigma.
# 9. Celithemis elisa (Calico Pennant): Belongs to the "pennants," known for being small perchers with weak, fluttery flight. This ecology makes a large pterostigma unnecessary.
# 10. Tholymis tillarga (Coral-tailed Cloudwing): A crepuscular species with a "dancing" flight. Research has confirmed it has reduced hindwing pterostigmata, likely for enhanced maneuverability.

# The indices for the species expected to have reduced pterostigmata are 9 and 10.
# The following code prints these indices as a comma-separated string.
species_with_reduced_pterostigmata = [9, 10]

# Print the final result in the specified format
output = ",".join(map(str, species_with_reduced_pterostigmata))
print(output)