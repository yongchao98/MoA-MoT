#
# Here is my plan:
# 1. Define the list of organisms and their life stages.
# 2. Identify which of these are mutualists of Asclepias fascicularis.
#    - Adult pollinators and defenders are mutualists.
#    - Larval herbivores are not mutualists.
#    - Larvae that do not interact with the plant are not mutualists.
# 3. Based on this, the mutualists are:
#    1) Danaus plexipus (Adult - pollinator)
#    2) Megachile frigidus (Adult - pollinator)
#    3) Formica rufa (Adult - nectar-feeder, potential defender)
#    4) Sphex ichneumoneus (Adult - pollinator)
#    5) Pepsis thisbe (Adult - pollinator)
#    6) Megachile ericetorum (Adult - pollinator)
# 4. Create a list containing the indices of these mutualists.
# 5. Print the indices separated by a comma.

mutualist_indices = [1, 2, 3, 4, 5, 6]

# Format the list of indices into a comma-separated string
output_string = ",".join(map(str, mutualist_indices))

# Print the final result
print(output_string)