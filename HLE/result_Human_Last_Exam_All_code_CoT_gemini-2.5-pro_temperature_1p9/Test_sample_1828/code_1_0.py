# The problem asks for the difference between the maximal and minimal possible
# cardinalities of the set X, where X is the set of cardinalities of uncountable
# maximal almost disjoint (MAD) families on omega.

# Step 1: Find the maximal possible cardinality of X.
# This occurs in a model of set theory where 2^omega = aleph_3,
# and where MAD families exist for all possible cardinalities between
# aleph_1 and aleph_3.
# The possible cardinalities are aleph_1, aleph_2, and aleph_3.
# So, in this model, X = {aleph_1, aleph_2, aleph_3}.
# The size of X is |X| = 3.
max_card_X = 3

# Step 2: Find the minimal possible cardinality of X.
# This occurs in a model of set theory where 2^omega = aleph_1 (the Continuum Hypothesis).
# In this model, all uncountable MAD families must have cardinality aleph_1.
# So, X = {aleph_1}.
# The size of X is |X| = 1.
min_card_X = 1

# Step 3: Calculate the difference.
difference = max_card_X - min_card_X

# Print the final calculation and result.
print(f"The maximal possible cardinality of X is {max_card_X}.")
print(f"The minimal possible cardinality of X is {min_card_X}.")
print(f"The difference is {max_card_X} - {min_card_X} = {difference}.")