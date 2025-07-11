# The task is to identify the indices of organisms that have a mutualistic relationship with Asclepias fascicularis.
# A mutualistic relationship is one where both organisms benefit.
# In this case, the plant benefits from pollination or defense, while the animal benefits from food (nectar).

# Step 1: Analyze each organism's relationship with the plant.
# - Adult pollinators (butterflies, bees, wasps) get nectar and pollinate the plant. This is mutualism.
# - Adult ants get nectar and can defend the plant from herbivores. This is mutualism.
# - Larvae are generally not mutualists in this context. Monarch larvae are herbivores, and the other larvae do not interact directly with the plant.

# Step 2: List the indices of the identified mutualists.
# - 1: Danaus plexipus (Adult) - Pollinator. Mutualist.
# - 2: Megachile frigidus (Adult) - Pollinator. Mutualist.
# - 3: Formica rufa (Adult) - Defender/Nectar-feeder. Mutualist.
# - 4: Sphex ichneumoneus (Adult) - Pollinator. Mutualist.
# - 5: Pepsis thisbe (Adult) - Pollinator. Mutualist.
# - 6: Megachile ericetorum (Adult) - Pollinator. Mutualist.
# - Larvae (7-12) are not mutualists.

# Step 3: Create a list of these indices.
mutualist_indices = [1, 2, 3, 4, 5, 6]

# Step 4: Format the list as a comma-separated string for the final output.
# We convert each number to a string and then join them with commas.
result_string = ",".join(map(str, mutualist_indices))

# Step 5: Print the result.
print(result_string)