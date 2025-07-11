import math

# Define the coordinates of the point of interest
x = 1
y = -1
z = 0

print(f"Calculating the magnetic field at point P(x,y,z) = ({x}, {y}, {z}).\n")

# --- Contribution from Wire 1 (along x-axis) ---
# The perpendicular distance from a point to the x-axis is sqrt(y^2 + z^2)
r1 = math.sqrt(y**2 + z**2)
# The magnitude of the magnetic field B1 is given by (mu_0 * I) / (2 * pi * r1).
# By the right-hand rule, its direction is -z.

print("--- Wire 1 (on x-axis) ---")
print(f"Perpendicular distance r1 = sqrt({y}^2 + {z}^2) = {r1}")
print(f"Magnitude |B1| = mu_0 * I / (2 * pi * {r1})")


# --- Contribution from Wire 2 (along y-axis) ---
# The perpendicular distance from a point to the y-axis is sqrt(x^2 + z^2)
r2 = math.sqrt(x**2 + z**2)
# The magnitude of the magnetic field B2 is given by (mu_0 * I) / (2 * pi * r2).
# By the right-hand rule, its direction is also -z.

print("\n--- Wire 2 (on y-axis) ---")
print(f"Perpendicular distance r2 = sqrt({x}^2 + {z}^2) = {r2}")
print(f"Magnitude |B2| = mu_0 * I / (2 * pi * {r2})")

# --- Total Magnetic Field ---
# Since both B1 and B2 are in the same direction, the total magnitude
# is the sum of the individual magnitudes.
# |B_total| = |B1| + |B2|
total_coefficient = (1/r1) + (1/r2)

print("\n--- Total Magnetic Field ---")
print("Since both fields point in the same direction (-z), we add their magnitudes.")
print("The final equation for the magnitude of the magnetic field is:")
print(f"|B_total| = |B1| + |B2| = (mu_0 * I / (2 * pi * {r1})) + (mu_0 * I / (2 * pi * {r2}))")
print(f"|B_total| = (mu_0 * I / (2 * pi)) * (1/{r1} + 1/{r2})")
print(f"|B_total| = (mu_0 * I / (2 * pi)) * ({total_coefficient})")
print(f"|B_total| = (mu_0 * I / pi)")
