import sympy
from sympy import sqrt, Symbol, pi

# Define 'd' as a symbolic variable to keep the answer in terms of 'd'.
d = Symbol('d', positive=True)

# --- Step 1: Determine the properties of the unrolled cone ---

# The base of the cone has a diameter 'd', so its radius 'r' is d/2.
r = d / 2

# The distance from the cone's apex to point P on the base circumference is 'd'.
# This is the cone's slant height 's'.
s = d

# When unrolled, the cone's surface becomes a sector of a circle.
# The radius of this sector is the slant height 's'.
R_sector = s

# The arc length of the sector is equal to the circumference of the cone's base.
L_arc = 2 * pi * r

# The angle 'theta' of the sector in radians is the arc length divided by the sector radius.
# We calculate it to understand the shape of the unrolled surface.
theta = L_arc / R_sector

print("--- Step 1: Unrolling the Cone ---")
print(f"The cone's base radius is r = d/2.")
print(f"The cone's slant height is s = {s}.")
print(f"The unrolled sector's radius is R = s = {s}.")
print(f"The unrolled sector's arc length is L = 2*pi*r = {sympy.simplify(L_arc)}.")
print(f"The sector's angle is theta = L/R = {sympy.simplify(theta)} radians, which is 180 degrees.")
print("Conclusion: The unrolled surface is a semicircle of radius d.")
print("-" * 20)

# --- Step 2: Finding the furthest point on the semicircle ---

print("\n--- Step 2: Locating the Furthest Point ---")
print("To visualize the problem, we place the semicircle in a 2D plane.")
print("The apex is at the origin (0, 0).")
print("The generator line passing through P is cut. This means P is represented by the two endpoints of the semicircle's diameter: (d, 0) and (-d, 0).")
print("The point 'Q' on the cone's surface that is furthest from 'P' must be equidistant from these two representations of P.")
print("This places Q on the line of symmetry of the semicircle, which is the y-axis (x=0).")
print("To maximize the distance from (d, 0) to a point Q=(0, y), we need to maximize y.")
print("The maximum value for y in the semicircle is d. So, the furthest point Q is at (0, d).")
print("-" * 20)

# --- Step 3: Calculating the maximum distance ---

# The coordinates of P can be taken as (d, 0) and Q as (0, d).
# We use the 2D distance formula: sqrt((x2-x1)^2 + (y2-y1)^2)
# This forms a right-angled triangle with both legs of length 'd'.

constant_factor = 2
print(f"\n--- Step 3: Final Calculation ---")
print("The distance is the hypotenuse of a right-angled triangle with legs of length d.")
print(f"Distance^2 = (d - 0)^2 + (0 - d)^2")
print(f"Distance^2 = d^2 + d^2 = {constant_factor} * d^2")
final_expression = f"d * sqrt({constant_factor})"
print(f"Distance = sqrt({constant_factor} * d^2) = {final_expression}")

print("\nTherefore, the furthest distance one can travel from point P is d times the square root of 2.")
print(f"<<<{final_expression}>>>")