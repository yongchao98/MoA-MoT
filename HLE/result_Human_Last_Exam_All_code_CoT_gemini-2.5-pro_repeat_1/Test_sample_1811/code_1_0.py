import math

# This script calculates the least number of zeros for a vector field on a
# compact manifold M with a non-empty boundary ∂M.

# The general formula is |χ(M) - (1/2)χ(∂M)|.

# Let's use a concrete example: the 2-dimensional disk, M = D².
# For the 2-disk, the Euler characteristic is χ(M) = 1.
# Its boundary is the circle, ∂M = S¹, for which the Euler characteristic is χ(∂M) = 0.

chi_M = 1
chi_dM = 0

# Calculate the result using the formula
# The result must be an integer. We use integer division // or convert to int.
result = int(abs(chi_M - chi_dM / 2))

# Print the general formula
print("The least number of zeros a vector field can have on a compact manifold M")
print("with a non-empty boundary ∂M is given by the formula:")
print("|χ(M) - (1/2) * χ(∂M)|")
print("\n--------------------------------------------------\n")

# Print the calculation for our example (M = D²)
print(f"For our example, the 2-disk M = D²:")
print(f"    χ(M) = {chi_M}")
print(f"    χ(∂M) = {chi_dM}")
print("\nPlugging these values into the formula:")
# The final equation with each number printed out
print(f"|{chi_M} - (1/2) * {chi_dM}| = |{chi_M} - {chi_dM/2}| = {result}")

print(f"\nSo, the least number of zeros a vector field can have on the 2-disk is {result}.")
