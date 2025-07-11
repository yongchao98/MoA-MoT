# This script identifies and prints the indices of species that perform photochemical synthesis.
# Based on the analysis, the following taxa utilize light to synthesize compounds:
# - 2, 6, 7: Classic oxygenic photosynthesis (algae, cyanobacteria, plant).
# - 9, 10, 12: Anoxygenic photosynthesis (bacteria).
# - 3: Vitamin D synthesis in humans.
# - 5: Light-driven ATP synthesis in archaea.

# The list of indices is determined from this biological knowledge.
photosynthetic_indices = [2, 3, 5, 6, 7, 9, 10, 12]

# The final output should be the indices separated by commas.
# The phrase "output each number in the final equation" is interpreted as
# printing the components of this final list.
output_string = ",".join(map(str, photosynthetic_indices))

print(output_string)