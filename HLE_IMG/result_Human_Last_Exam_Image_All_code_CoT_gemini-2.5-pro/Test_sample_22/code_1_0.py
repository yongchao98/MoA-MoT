# Based on knot theory analysis, we determine for each knot if it's an unknot.

# K_1: Not the unknot (it's the 6_1 knot).
# K_2: The unknot (a known complex diagram).
# K_3: The unknot (simplifies with one Reidemeister II move).
# K_4: Not the unknot (it's the figure-eight knot).
# K_5: The unknot (simplifies with a sequence of Reidemeister moves).
# K_6: Not the unknot (it's the 10_161 knot).

# We represent these findings as a list of boolean values.
is_unknot_list = [False, True, True, False, True, False]

# We will now generate the list of 1-based indices for the unknots.
unknot_indices = []

# Iterate through the boolean list with indices starting from 1.
for i, is_unknot in enumerate(is_unknot_list):
    # The knot index is i + 1
    knot_index = i + 1
    if is_unknot:
        unknot_indices.append(knot_index)

# Print the final list of indices.
print(unknot_indices)