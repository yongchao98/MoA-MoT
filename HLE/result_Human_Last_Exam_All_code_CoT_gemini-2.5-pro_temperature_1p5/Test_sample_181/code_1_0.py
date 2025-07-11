import math

# The problem is to calculate foo(7). Our analysis shows that foo(n) counts the
# number of "right group" structures on a labeled set of n elements.
# For a prime number n, a right group must be either a group of order n
# or a right-zero semigroup of order n.

n = 7

# Case 1: The structure is a group of order 7.
# Since 7 is a prime number, any group of order 7 is isomorphic to the cyclic group Z_7.
# The number of distinct group structures (i.e., distinct multiplication tables) on a
# labeled set of 7 elements is given by (7-1)!.
num_groups = math.factorial(n - 1)

# Case 2: The structure is a right-zero semigroup of order 7.
# A right-zero semigroup has the operation defined as x * y = y. For a labeled set,
# this operation corresponds to exactly one unique multiplication table.
# For example, for elements {1, 2, ..., 7}, the (i, j)-th entry of the table is always j.
num_rz_semigroups = 1

# The total number of structures is the sum of the counts from these two disjoint cases.
total = num_groups + num_rz_semigroups

# As requested, printing each number in the final equation.
print(f"{num_groups} + {num_rz_semigroups} = {total}")