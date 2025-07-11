import math

# The problem asks for the magnitude of the magnetic field at the point (x,y,z)=(1,-1,0).

# Let's define the coordinates of the point P.
x, y, z = 1, -1, 0

# --- Step 1: Analyze the contribution from Wire 1 (on the x-axis) ---

# The perpendicular distance 'r1' from a point (x, y, z) to the x-axis is sqrt(y^2 + z^2).
r1 = math.sqrt(y**2 + z**2)

# By the right-hand rule, with current in the +x direction, the magnetic field at a
# point with a negative y-coordinate (like P) points in the positive z-direction.
# So, the vector B1 is B1 = (mu_0*I / (2*pi*r1)) * k_hat, where k_hat is the unit vector in the z-direction.


# --- Step 2: Analyze the contribution from Wire 2 (on the y-axis) ---

# The perpendicular distance 'r2' from a point (x, y, z) to the y-axis is sqrt(x^2 + z^2).
r2 = math.sqrt(x**2 + z**2)

# By the right-hand rule, with current in the +y direction, the magnetic field at a
# point with a positive x-coordinate (like P) points in the negative z-direction.
# So, the vector B2 is B2 = - (mu_0*I / (2*pi*r2)) * k_hat.


# --- Step 3: Calculate the total field and its magnitude ---

# The total magnetic field is the vector sum B_total = B1 + B2.
# Since both vectors are along the z-axis, we can sum their components:
# B_total_z = (mu_0*I / (2*pi*r1)) - (mu_0*I / (2*pi*r2))
# B_total_z = (mu_0*I / (2*pi)) * (1/r1 - 1/r2)
# The magnitude is the absolute value of this component.

# Let's print the final equation with the numbers plugged in.
# The term (μ₀·I / 2π) is a common constant factor.
constant_factor_str = "(\u03BC\u2080\u00B7I / 2\u03C0)"

print(f"The magnitude of the total magnetic field |B_total| is given by the superposition of the fields from each wire.")
print(f"|B_total| = | B_wire1_z + B_wire2_z |")
print(f"|B_total| = | {constant_factor_str}/r1 - {constant_factor_str}/r2 |")
print("\nFirst, we calculate the distances r1 and r2 for the point (1, -1, 0):")
print(f"r1 = \u221A(y\u00B2 + z\u00B2) = \u221A(({y})\u00B2 + {z}\u00B2) = {r1}")
print(f"r2 = \u221A(x\u00B2 + z\u00B2) = \u221A({x}\u00B2 + {z}\u00B2) = {r2}")

print("\nSubstituting these distances into the final equation for the magnitude:")
# The prompt requires showing each number in the final equation.
# Here we show the contribution from each term explicitly.
term1_coeff = 1 / r1
term2_coeff = 1 / r2
final_coefficient = term1_coeff - term2_coeff

print(f"|B_total| = {constant_factor_str} * | (1/r1) - (1/r2) |")
print(f"|B_total| = {constant_factor_str} * | (1/{r1}) - (1/{r2}) |")
print(f"|B_total| = {constant_factor_str} * | {term1_coeff} - {term2_coeff} |")
print(f"|B_total| = {constant_factor_str} * |{final_coefficient}|")
print(f"|B_total| = 0")