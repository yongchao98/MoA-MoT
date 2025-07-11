# Based on the analysis of the life habits of the immatures of each tribe,
# this script identifies those unlikely to be collected with a beat-sheet method.

# Unlikely tribes are those whose immatures do not live freely on plant foliage.
# 1. Apis: Larvae live in hives.
# 3. Eupholini: Larvae are internal plant borers.
# 4. Acritini: Larvae live in decaying matter.
# 7. Acanthocerini: Larvae live in decaying wood/fungus.

# Create a list of the indices for these unlikely tribes.
unlikely_indices = [1, 3, 4, 7]

# The indices are already in ascending order.
# The following code will format and print the result as a comma-separated string.
print(','.join(map(str, unlikely_indices)))