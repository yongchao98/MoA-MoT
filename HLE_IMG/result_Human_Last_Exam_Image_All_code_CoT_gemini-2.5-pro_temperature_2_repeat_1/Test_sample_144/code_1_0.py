import collections

# Step 1: Define the identified mimics and damage-causers.
# This is based on the visual analysis of the six panels.
mimics = ["A", "C", "E"]
causers = ["B", "D", "F"]

# Step 2: Define the pairings based on the logic described.
# Each key is a mimic, and its value is the corresponding damage-causer.
# The reasoning for each pair is:
# A -> B: The stripe on the adult beetle (A) mimics the linear feeding damage of the larva (B).
# C -> F: The tattered, decayed look of the moth's wings (C) mimics the ragged damage caused by the chewing katydid (F).
# E -> D: The veiny appearance of the leaf insect (E) mimics a skeletonized leaf, a type of damage caused by the leaf beetle (D).
pairings = {
    "A": "B",
    "C": "F",
    "E": "D"
}

# Step 3: Create the output string in the format "AB, CD, EF".
# We iterate through a defined order to ensure consistency, though any order of pairs is valid.
output_pairs = []
# We use a defined order for printing consistency.
mimic_order = ["A", "C", "E"]
for mimic in mimic_order:
  causer = pairings[mimic]
  # For each pair, the output string like "AB" is formed.
  # The explanation for each part of the 'equation' is:
  # {mimic} represents the insect that mimics leaf damage.
  # {causer} represents the insect that causes the leaf damage being mimicked.
  output_pairs.append(f"{mimic}{causer}")

# Step 4: Print the final result.
print(", ".join(output_pairs))