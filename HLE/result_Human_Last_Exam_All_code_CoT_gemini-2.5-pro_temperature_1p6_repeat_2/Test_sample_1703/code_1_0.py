import math

# This script derives and prints the formula for the area of the triangle T(t).

# --- Problem Parameters ---
# Radius of the circle C, which is equal to the side length of the inscribed hexagon H.
hexagon_side_length = 10
# Speed of the triangle's vertices along the hexagon's sides.
vertex_speed = 1

# --- Derivation of Formula Components ---

# At t=0, the triangle's vertices are at the midpoints of alternating sides.
# The resulting triangle T(0) is equilateral, and its side length is 1.5 times the hexagon's side length.
initial_triangle_side_length = 1.5 * hexagon_side_length

# The constant term in the area formula comes from the initial area.
# It's related to the square of the initial side length.
# This term is (1.5 * 10)^2 = 15^2
constant_term = initial_triangle_side_length**2

# The time-dependent term arises from the movement of the vertices.
# The squared side length of the triangle T(t) is s_T(t)^2 = (15)^2 + 3 * (v*t)^2
# Since v=1, the coefficient for the t^2 term is 3.
time_dependent_coefficient = 3

# The area of an equilateral triangle is (sqrt(3)/4) * (side_length)^2.
# These are the components from this standard formula.
sqrt_argument = 3
denominator = 4
time_power = 2

# --- Final Formula Output ---
# The prompt asks to output each number in the final equation.
# The final formula is A(t) = (sqrt(3)/4) * (225 + 3*t^2)

print("The area of the triangle T(t) as a function of time t is given by the formula:")
print(f"A(t) = (sqrt({sqrt_argument}) / {denominator}) * ({int(constant_term)} + {time_dependent_coefficient}*t**{time_power})")

# For clarity, here is the simplified expression with approximate numerical values.
coeff_constant = (math.sqrt(sqrt_argument) / denominator) * constant_term
coeff_t_squared = (math.sqrt(sqrt_argument) / denominator) * time_dependent_coefficient
print("\nIn a simplified numerical form, the formula is approximately:")
print(f"A(t) ≈ {coeff_constant:.4f} + {coeff_t_squared:.4f} * t**{time_power}")
