import math

# Let M be a compact manifold with a non-empty boundary ∂M.
# χ(M) is the Euler characteristic of M.
# χ(∂M) is the Euler characteristic of its boundary ∂M.

# The user should replace these placeholder values with the Euler
# characteristics of the specific manifold they are interested in.
# As an example, we use a 2-dimensional disk (D^2).
# For a disk, M = D^2, the Euler characteristic is χ(M) = 1.
# Its boundary is a circle, ∂M = S^1, for which χ(∂M) = 0.
chi_M = 1
chi_dM = 0

# The Poincaré-Hopf theorem, generalized for manifolds with boundary,
# states that the sum of the indices of the zeros of a vector field is
# a topological invariant given by the formula: χ(M) - χ(∂M)/2.
#
# A key property is that for any compact manifold M, the Euler
# characteristic of its boundary, χ(∂M), is always an even number.
# Therefore, χ(∂M)/2 is always an integer.

# We use integer division // as χ(∂M) is always even.
half_chi_dM = chi_dM // 2

# Calculate the sum of the indices
index_sum = chi_M - half_chi_dM

# The least number of zeros a vector field can have is the
# absolute value of this sum.
min_zeros = abs(index_sum)

# Output the formula and the calculation.
print("The formula for the minimum number of zeros of a vector field on a manifold M with boundary ∂M is |χ(M) - χ(∂M)/2|.")
print(f"For the example case where χ(M) = {chi_M} and χ(∂M) = {chi_dM}:")
print(f"Least number of zeros = |{chi_M} - {chi_dM}/2| = |{chi_M} - {half_chi_dM}| = |{index_sum}| = {min_zeros}")

# Let's try another example: a 3-dimensional ball (B^3).
# For M = B^3, χ(M) = 1.
# Its boundary is a 2-sphere, ∂M = S^2, for which χ(∂M) = 2.
chi_M_ball = 1
chi_dM_ball = 2
half_chi_dM_ball = chi_dM_ball // 2
index_sum_ball = chi_M_ball - half_chi_dM_ball
min_zeros_ball = abs(index_sum_ball)

print("\nFor the example case of a 3-ball where χ(M) = {chi_M_ball} and χ(∂M) = {chi_dM_ball}:")
print(f"Least number of zeros = |{chi_M_ball} - {chi_dM_ball}/2| = |{chi_M_ball} - {half_chi_dM_ball}| = |{index_sum_ball}| = {min_zeros_ball}")
