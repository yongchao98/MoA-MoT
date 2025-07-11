# This script calculates the difference between the maximal and minimal
# possible cardinalities of the set X, based on set-theoretic results.

# Step 1: Define the maximal possible cardinality of X.
# From the analysis, under the given hypotheses (not CH and 2^{omega_1}=omega_3),
# the set of possible cardinalities for MAD families, X, is a subset of
# {omega_1, omega_2, omega_3}.
# It is consistent with ZFC to construct a model where MAD families of all these
# three cardinalities exist. In such a model, the size of X is 3.
max_card_X = 3

# Step 2: Define the minimal possible cardinality of X.
# The size of X must be at least 1.
# It is consistent with ZFC to construct a model where all MAD families
# have the same cardinality (e.g., omega_2). In such a model, the size of X is 1.
min_card_X = 1

# Step 3: Calculate the difference.
difference = max_card_X - min_card_X

print(f"The maximal possible cardinality of the set X is: {max_card_X}")
print(f"The minimal possible cardinality of the set X is: {min_card_X}")
print(f"The difference is the maximal value minus the minimal value.")
print(f"Difference = {max_card_X} - {min_card_X} = {difference}")
