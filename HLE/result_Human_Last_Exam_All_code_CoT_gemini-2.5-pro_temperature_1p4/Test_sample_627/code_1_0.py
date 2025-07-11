# Program to calculate the braid index upper bound for the three-twist knot using Vogel's algorithm.

# 1. Define the properties of the three-twist knot (6_1).

# c: The minimal crossing number for the three-twist knot.
c = 6

# For any knot diagram, the total number of regions created on the plane is c + 2.
num_total_regions = c + 2

# n: The number of regions of a single color in a checkerboard coloring.
# This is half of the total regions.
n = num_total_regions // 2

# k: The number of connected components of the checkerboard graph.
# For the minimal diagram of an alternating knot like the three-twist knot,
# the checkerboard graph is connected.
k = 1

# 2. Apply Vogel's algorithm formula to calculate the upper bound.
# Formula: Upper Bound = c - n + k + 1
upper_bound = c - n + k + 1

# 3. Print the explanation and the result.
print("This script calculates an upper bound for the braid index of the three-twist knot (6_1) using Vogel's algorithm.")
print("The formula for the upper bound is: c - n + k + 1\n")
print("For the three-twist knot:")
print(f"  - The crossing number is c = {c}.")
print(f"  - The number of regions of a single color is n = {n}.")
print(f"  - The number of connected components of the checkerboard graph is k = {k}.\n")

print("Calculating the upper bound:")
# The final output needs to show each number in the equation.
print(f"Upper Bound = {c} - {n} + {k} + 1 = {upper_bound}")