# The task is to identify the mutualists of Asclepias fascicularis from a given list.
# A mutualistic relationship is one where both species benefit.
# For plants and insects, this usually involves pollination or protection in exchange for nectar.

# Analysis:
# 1) Danaus plexipus (Adult Monarch Butterfly): Pollinator. Feeds on nectar. Mutualist.
# 2) Megachile frigidus (Adult Leafcutter Bee): Pollinator. Feeds on nectar and pollen. Mutualist.
# 3) Formica rufa (Adult Red Wood Ant): Nectar feeder. Can protect the plant from herbivores. Mutualist.
# 4) Sphex ichneumoneus (Adult Great Golden Digger Wasp): Pollinator. Feeds on nectar. Mutualist.
# 5) Pepsis thisbe (Adult Tarantula Hawk Wasp): Pollinator. Feeds on nectar. Mutualist.
# 6) Megachile ericetorum (Adult Leafcutter Bee): Pollinator. Feeds on nectar and pollen. Mutualist.
# 7) Danaus plexipus (Larva/Caterpillar): Herbivore. Eats milkweed leaves. Not a mutualist.
# 8) Megachile frigidus (Larva): Non-interactive. Develops in a nest. Not a mutualist.
# 9) Formica rufa (Larva): Non-interactive. Develops in a nest. Not a mutualist.
# 10) Sphex ichneumoneus (Larva): Non-interactive. Develops on prey in a burrow. Not a mutualist.
# 11) Pepsis thisbe (Larva): Non-interactive. Develops on a paralyzed tarantula. Not a mutualist.
# 12) Megachile ericetorum (Larva): Non-interactive. Develops in a nest. Not a mutualist.

# Collect the indices of all the identified mutualists.
mutualist_indices = [1, 2, 3, 4, 5, 6]

# Prepare the output string.
# We convert each integer index to a string and then join them with a comma.
output_string = ",".join(map(str, mutualist_indices))

# Print the final result.
print(output_string)