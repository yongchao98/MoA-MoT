# The problem is to find the minimum number of guards to observe the exterior
# of a unit ball in n-dimensional space, with guards on the boundary.

# This problem is a known case of the exterior illumination problem.
# A theorem by Zalgaller states that for a convex body of constant width in R^n,
# n+1 guards on the boundary are sufficient.
# A sphere is a body of constant width.

# We are in 3D space, so the dimension 'n' is 3.
n = 3

# The number of guards required is n + 1.
num_guards = n + 1

# Print the final equation and the result.
# The format is number = number + number.
print(f"The dimension of the space is n = {n}.")
print("The minimum number of guards required is n + 1.")
print("So, the final calculation is:")
print(f"{num_guards} = {n} + 1")
