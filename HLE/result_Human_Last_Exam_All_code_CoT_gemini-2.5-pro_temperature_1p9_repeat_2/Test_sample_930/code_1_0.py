# The task is to identify the mutualists of Asclepias fascicularis from a given list.
# A mutualistic relationship benefits both species. In this context, it primarily refers to pollination.

# 1. Adult Danaus plexipus (Monarch Butterfly) is a pollinator. Mutualist.
# 2. Adult Megachile frigida (Leaf-cutter Bee) is a pollinator. Mutualist.
# 3. Adult Formica rufa (Ant) is a nectar robber, not a reliable pollinator. Not a mutualist.
# 4. Adult Sphex ichneumoneus (Digger Wasp) is a pollinator. Mutualist.
# 5. Adult Pepsis thisbe (Tarantula Hawk Wasp) is a pollinator. Mutualist.
# 6. Adult Megachile ericetorum (Leaf-cutter Bee) is a pollinator. Mutualist.
# 7. Larval Danaus plexipus is a herbivore, which harms the plant. Not a mutualist.
# 8-12. The larvae of bees, ants, and wasps do not interact directly with the plant. Not mutualists.

# Therefore, the indices of the mutualists are 1, 2, 4, 5, and 6.
mutualist_indices = [1, 2, 4, 5, 6]

# Convert the list of numbers to a comma-separated string for printing.
# Using map(str, ...) to convert each integer to a string before joining.
final_answer = ",".join(map(str, mutualist_indices))

print(final_answer)