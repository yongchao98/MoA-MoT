# The plan is to identify the organisms that have a mutualistic relationship with Asclepias fascicularis.
# Mutualists benefit the plant, primarily through pollination or defense, while also receiving a benefit, like nectar.
# 1) Danaus plexipus (Adult): Pollinator. Mutualist.
# 2) Megachile frigidus (Adult): Pollinator. Mutualist.
# 3) Formica rufa (Adult): Nectar feeder and potential defender. Mutualist.
# 4) Sphex ichneumoneus (Adult): Pollinator. Mutualist.
# 5) Pepsis thisbe (Adult): Pollinator. Mutualist.
# 6) Megachile ericetorum (Adult): Pollinator. Mutualist.
# 7-12) All larvae listed are not mutualists. They are either herbivores (Danaus plexipus) or do not interact directly with the plant.

# Create a list of the indices for the mutualists.
mutualist_indices = [1, 2, 3, 4, 5, 6]

# The output should be a comma-separated string of these indices.
# We convert each number to a string and then join them with a comma.
output_string = ",".join(map(str, mutualist_indices))

# Print the final result.
print(output_string)