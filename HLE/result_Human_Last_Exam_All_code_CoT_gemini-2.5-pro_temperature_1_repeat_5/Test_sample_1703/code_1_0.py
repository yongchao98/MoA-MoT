import math

# Plan:
# 1. Recognize that the hexagon's rotation is irrelevant to the triangle's area.
# 2. Determine that the triangle T(t) is always equilateral due to symmetry.
# 3. Derive the area formula in terms of the position of the vertices on the hexagon's sides.
# 4. Express the vertex position as a function of time t.
# 5. Combine these to get the final area formula Area(t).
# 6. The code will print the components of this final formula as requested.

# Based on the derivation, the area of the triangle T(t) can be expressed as:
# Area(t) = A + B * [u(t)]^2
# where u(t) is a periodic function representing the displacement of the vertices
# from the midpoint of their respective hexagon sides.

# The values for A and B are derived from the geometry of the problem (R=10).
# A = (225 * sqrt(3)) / 4
# B = (3 * sqrt(3)) / 4

# The function u(t) is a triangular wave because the vertices move with constant speed
# back and forth along the hexagon's sides (length 10).
# The movement from midpoint to one end (5 units) takes 5s. The full cycle
# (midpoint -> end -> other end -> midpoint) takes 20s. The amplitude of displacement
# from the midpoint is 5.

# Define the numerical components of the final equation
A_numerator = 225
A_denominator = 4
A_value = (A_numerator / A_denominator) * math.sqrt(3)

B_numerator = 3
B_denominator = 4
B_value = (B_numerator / B_denominator) * math.sqrt(3)

period = 20
amplitude = 5

print("The area of the triangle T(t) is given by the function:")
print(f"Area(t) = A + B * [u(t)]^2")
print("\nWhere the numbers in the final equation are:")
print(f"Constant term A = ({A_numerator} * sqrt(3)) / {A_denominator}")
print(f"                 (which is approximately {A_value:.4f})")
print(f"Coefficient B = ({B_numerator} * sqrt(3)) / {B_denominator}")
print(f"                (which is approximately {B_value:.4f})")
print("\nAnd u(t) is a periodic triangular wave function describing the vertex displacement from the side's midpoint.")
print(f"The function u(t) has a period of {period} seconds and an amplitude of {amplitude} units.")
