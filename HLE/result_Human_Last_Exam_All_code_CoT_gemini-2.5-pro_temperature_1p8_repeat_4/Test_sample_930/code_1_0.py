# The plan is to identify the organisms that have a mutualistic relationship
# with Asclepias fascicularis, meaning both the plant and the animal benefit.
# In this scenario, the mutualists are the adult insects that pollinate the plant
# while feeding on its nectar.

# Based on ecological roles:
# - Adult butterflies, bees, and wasps that visit flowers for nectar are pollinators.
# - The Danaus plexippus larva is an herbivore, damaging the plant.
# - Other larvae (bees, ants, parasitic wasps) do not interact with the plant directly.
# - Ants often rob nectar or protect herbivores (aphids), which is not mutualistic.

# Therefore, the mutualists are the adult forms of Danaus plexipus (1),
# Megachile frigidus (2), Sphex ichneumoneus (4), Pepsis thisbe (5),
# and Megachile ericetorum (6).

mutualist_indices = [1, 2, 4, 5, 6]

# Convert the list of numbers to a comma-separated string for the output.
output = ",".join(map(str, mutualist_indices))

print(output)