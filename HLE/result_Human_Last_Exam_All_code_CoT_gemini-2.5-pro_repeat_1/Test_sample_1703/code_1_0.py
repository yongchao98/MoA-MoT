import math

# This script calculates the area of the triangle T(t) as a function of time t.

# --- Problem Parameters ---
# Radius of the circumscribing circle for the hexagon
R = 10
# Speed of the triangle's vertices along the hexagon's sides
v = 1

# --- Geometric Derivation ---

# 1. The side length of the hexagon is equal to the radius R.
# The distance from the center of the hexagon to the midpoint of a side (apothem, h)
# is calculated as h = R * sqrt(3) / 2.
# We will use the square of this distance: h^2 = (3/4) * R^2.
h_squared = (3/4) * (R**2)

# 2. At time t, a vertex of the triangle has moved a distance d = v*t from the
# midpoint of a hexagon side. Since v=1, d=t. The square of this distance is d^2 = t^2.

# 3. The distance of a triangle vertex from the center of the hexagon (d_vertex) can be
# found using the Pythagorean theorem: d_vertex^2 = h^2 + d^2.
# d_vertex^2(t) = (3/4)*R^2 + t^2

# 4. For an equilateral triangle with side length s, the distance from its center
# to a vertex is d_vertex = s / sqrt(3). Thus, the side length squared is s^2 = 3 * d_vertex^2.
# s(t)^2 = 3 * ((3/4)*R^2 + t^2) = (9/4)*R^2 + 3*t^2.
constant_part_s_squared = (9/4) * (R**2)
t_squared_coeff_s_squared = 3 * (v**2)

# 5. The area of an equilateral triangle is Area = (sqrt(3)/4) * s^2.
# Substituting our expression for s(t)^2:
# Area(t) = (sqrt(3)/4) * ((9/4)*R^2 + 3*t^2)
# We can now substitute the numerical value R=10.
final_constant_term = int(constant_part_s_squared)
final_t_coeff = int(t_squared_coeff_s_squared)

# --- Final Result ---
print("The area of the triangle T(t), denoted as A(t), is derived from its side length.")
print("The key steps are realizing the triangle remains equilateral and its area is independent of the hexagon's rotation.")
print("\nThe final formula for the area as a function of time t is:")
print(f"A(t) = (sqrt(3) / 4) * ({final_constant_term} + {final_t_coeff} * t^2)")