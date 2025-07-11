# Step 1: Analyze the list to identify mutualists of Asclepias fascicularis.
# A mutualist relationship is one where both species benefit.

# For Asclepias fascicularis (Narrowleaf Milkweed), mutualists are typically pollinators or defenders.

# Adults:
# 1) Danaus plexipus (Monarch butterfly): Adult butterflies are pollinators. They drink nectar and transfer pollen. This is a mutualistic relationship.
# 2) Megachile frigidus (Leaf-cutter bee): Adult bees are pollinators. This is a mutualistic relationship.
# 3) Formica rufa (Red wood ant): Ants are often attracted to extrafloral nectaries on milkweeds, and in return for the food, they defend the plant from herbivores. This is a mutualistic relationship.
# 4) Sphex ichneumoneus (Great golden digger wasp): Adult wasps drink nectar and act as pollinators. This is a mutualistic relationship.
# 5) Pepsis thisbe (Tarantula hawk wasp): Adult wasps are important pollinators for many plants, including milkweeds. This is a mutualistic relationship.
# 6) Megachile ericetorum (Leaf-cutter bee): Like other bees, adults are pollinators. This is a mutualistic relationship.

# Larvae:
# 7) Danaus plexipus (Monarch caterpillar): The larva is a herbivore that eats the milkweed's leaves. This is a harmful relationship for the plant, not mutualistic.
# 8) Megachile frigidus (Bee larva): The larva develops inside a nest and is fed by the adult; it does not interact with the plant. Not mutualistic.
# 9) Formica rufa (Ant larva): The larva develops inside the nest; it does not interact with the plant. Not mutualistic.
# 10) Sphex ichneumoneus (Wasp larva): The larva is a parasitoid, feeding on another insect provided by the adult. It does not interact with the plant. Not mutualistic.
# 11) Pepsis thisbe (Wasp larva): The larva is a parasitoid of tarantulas. It does not interact with the plant. Not mutualistic.
# 12) Megachile ericetorum (Bee larva): The larva develops inside a nest; it does not interact with the plant. Not mutualistic.

# Step 2: Compile the indices of the mutualists.
mutualist_indices = [1, 2, 3, 4, 5, 6]

# Step 3: Format the output as a comma-separated string.
# Convert each integer in the list to a string
string_indices = [str(index) for index in mutualist_indices]

# Join the list of strings with a comma
output_string = ",".join(string_indices)

# Step 4: Print the final result.
print(output_string)