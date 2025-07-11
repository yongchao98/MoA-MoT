# The task is to identify dragonfly species from a list that are expected to have reduced pterostigmata based on their ecology.
# A reduced pterostigma is strongly associated with a highly migratory, gliding lifestyle, as seen in pelagic (open-air) species.
# We will evaluate the ecology of each species to identify the long-distance migrants.

# 1) Didymops transversa: Stream patroller, not a long-distance migrant.
# 2) Urothemis edwarsi: A percher, not a long-distance migrant.
# 3) Macrodiplax balteata: A known coastal and long-distance migrant. This ecology suggests reduced pterostigmata.
# 4) Pantala flavescens: The "Globe Skimmer," a famous trans-oceanic migrant that spends most of its life gliding. This is a classic example of a species with reduced pterostigmata.
# 5) Orthetrum cancellatum: A percher, not a long-distance migrant.
# 6) Libelulla quadrimaculata: Can be migratory, but is a strong flier with prominent pterostigmata, not primarily a glider.
# 7) Libelulla pulchela: A percher, not a long-distance migrant.
# 8) Sympetrum corruptum: A migrant, but not as specialized for gliding as others; its pterostigmata are not notably reduced.
# 9) Celithemis elisa: A percher, not a migrant.
# 10) Tholymis tillarga: A crepuscular, migratory species known for its gliding flight. This ecology suggests reduced pterostigmata.

# Based on this analysis, the species with the appropriate ecologies are 3, 4, and 10.
species_indices_with_reduced_pterostigmata = [3, 4, 10]

# Convert the list of numbers to a list of strings for joining
indices_as_strings = [str(index) for index in species_indices_with_reduced_pterostigmata]

# Join the strings with a comma to create the final output string
final_answer = ",".join(indices_as_strings)

print(final_answer)