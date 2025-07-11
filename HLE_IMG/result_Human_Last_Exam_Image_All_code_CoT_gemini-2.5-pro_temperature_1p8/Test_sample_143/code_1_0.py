import collections

# Step 1: Analyze each image pair to determine the status of the species.
# Pair A shows Garlic Mustard. The right image is a dense monoculture, a hallmark of invasion. The left is a mixed community, typical of a native range.
# This corresponds to option 3: 'right invaded'.
analysis_A = 3

# Pair B shows Bigleaf Lupine. The left image shows a dense, dominant stand characteristic of its invasive behavior in places like Iceland. The right is much sparser.
# This corresponds to option 4: 'left invaded'.
analysis_B = 4

# Pair C shows Papaya. It's a species native to the Americas but widely naturalized in tropics worldwide. It's plausible both photos are from outside its native range where it has become established.
# This corresponds to option 2: 'both invasive'.
analysis_C = 2

# Step 2: Assemble the final answer.
# The final answer is the sequence of indices for A, B, and C.
final_answer = f"{analysis_A}, {analysis_B}, {analysis_C}"

# Step 3: Print the final answer.
print(final_answer)