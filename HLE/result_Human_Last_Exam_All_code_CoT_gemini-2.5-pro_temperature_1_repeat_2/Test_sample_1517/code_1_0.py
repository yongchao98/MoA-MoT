import math

# For any integer n, it is possible to construct a countably infinite 
# non-discrete Hausdorff topological group G such that the minimum index 
# I_G is equal to n.
#
# Let n be an arbitrarily large integer.
# We can find a group G_n where I_{G_n} = n.
# For example, if we let n = 1000, there is a group G_1000 with I_{G_1000} = 1000.
# If we let n = 1000000, there is a group G_1000000 with I_{G_1000000} = 1000000.
#
# This means the set of possible values for I_G is {1, 2, 3, ..., n, ...} and also includes infinite values like aleph_0.
# Since this set is not bounded by any finite number, there is no largest finite value.
# The supremum of this set of values is infinity.

largest_value = math.inf
print(largest_value)