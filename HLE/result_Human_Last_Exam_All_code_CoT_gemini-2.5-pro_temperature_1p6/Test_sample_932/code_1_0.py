# The task is to identify which of the given insect tribes' immatures
# are unlikely to be collected using a beat-sheet method.

# Based on entomological knowledge:
# 1) Apis: Immatures are in a hive, not on foliage. -> Unlikely.
# 2) Melipotini: Caterpillars feed on foliage. -> Likely.
# 3) Eupholini: Weevil larvae are internal borers. -> Unlikely.
# 4) Acritini: Larvae live in decaying matter, not on foliage. -> Unlikely.
# 5) Oxyptilini: Caterpillars feed on foliage/flowers. -> Likely.
# 6) Dictyophorini: Nymphs live on plant surfaces. -> Likely.
# 7) Acanthocerini: Larvae live in rotting wood. -> Unlikely.

# The indices of the unlikely tribes are 1, 3, 4, and 7.
unlikely_indices = [1, 3, 4, 7]

# Sort the indices in ascending order, though they already are.
unlikely_indices.sort()

# Convert the list of numbers into a comma-separated string for the output.
result = ",".join(map(str, unlikely_indices))

# Print the final result.
print(result)