import math

# An icosahedron has 12 vertices. The axes of 5-fold rotational symmetry
# pass through pairs of opposite vertices.
num_vertices = 12
num_axes = num_vertices // 2

# For each axis of 5-fold symmetry, there are φ(5) rotations that have order 5.
# φ(n) is Euler's totient function, which counts the positive integers up to n
# that are relatively prime to n. math.phi is available from Python 3.8 onwards.
# For older versions, phi(5) = 5 * (1 - 1/5) = 4.
num_rotations_per_axis = 4 # math.phi(5)

# The total number of elements of order 5 in the rotational group A_5 is
# the number of axes multiplied by the number of order-5 rotations per axis.
# This corresponds to the number of elements in H3 that satisfy the given conditions.
total_elements = num_axes * num_rotations_per_axis

print("The number of elements is the result of the following equation:")
print(f"Number of axes through opposite vertices * Number of order-5 rotations per axis = Total elements")
print(f"{num_axes} * {num_rotations_per_axis} = {total_elements}")
