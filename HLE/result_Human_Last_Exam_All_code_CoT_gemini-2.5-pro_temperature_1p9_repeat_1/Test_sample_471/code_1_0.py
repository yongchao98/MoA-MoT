# The minimal number of critical points for a smooth function on a compact manifold
# is bounded below by the sum of its Betti numbers.

# Betti numbers for the 2-torus (T^2):
b_0 = 1  # number of connected components
b_1 = 2  # number of 1-dimensional ("circular") holes
b_2 = 1  # number of 2-dimensional ("void") holes

# Calculate the sum of the Betti numbers.
# This sum gives the minimal number of critical points.
minimal_critical_points = b_0 + b_1 + b_2

print("The minimal number of critical points is the sum of the Betti numbers of the 2-torus.")
print("The Betti numbers are b_0, b_1, and b_2.")
print(f"b_0 = {b_0}")
print(f"b_1 = {b_1}")
print(f"b_2 = {b_2}")
print("\nThe final calculation is:")
# Output each number in the final equation
print(f"{b_0} + {b_1} + {b_2} = {minimal_critical_points}")