# Step 1: Determine the maximal possible cardinality of X.
# Based on set-theoretic consistency results, it is possible to have a model of ZFC
# satisfying the problem's conditions where the set of cardinalities of uncountable
# maximal almost disjoint families is X = {omega_1, omega_2}.
# In this model, the cardinality of X is 2. This is the maximum possible value.
max_possible_cardinality_of_X = 2

# Step 2: Determine the minimal possible cardinality of X.
# It is also consistent with ZFC to have a model where the only possible size for an
# uncountable maximal almost disjoint family is omega_2 (the value of the continuum).
# In this model, the set of such cardinalities is X = {omega_2}.
# The cardinality of X is 1. This is the minimum possible value, as X cannot be empty.
min_possible_cardinality_of_X = 1

# Step 3: Calculate the difference.
difference = max_possible_cardinality_of_X - min_possible_cardinality_of_X

# Step 4: Output the result in the specified format.
print(f"The maximal possible cardinality of X is {max_possible_cardinality_of_X}.")
print(f"The minimal possible cardinality of X is {min_possible_cardinality_of_X}.")
print("The difference is:")
print(f"{max_possible_cardinality_of_X} - {min_possible_cardinality_of_X} = {difference}")