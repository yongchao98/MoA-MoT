# The species with ecologies based on extensive gliding and long-distance migration are expected to have reduced pterostigmata.
# This adaptation is most prominent in the "glider guild" of dragonflies.
# - Pantala flavescens (4) is known as the "Wandering Glider" and is the most famous example.
# - Tholymis tillarga (10) is another well-known migratory glider.
# These two species are the strongest candidates.

# Define the list of indices for the identified species.
species_with_reduced_pterostigmata = [4, 10]

# Format the output as a comma-separated string.
# The map() function converts each integer in the list to a string.
# The join() method concatenates the string elements with a comma.
output_string = ",".join(map(str, species_with_reduced_pterostigmata))

# Print the final result.
print(output_string)