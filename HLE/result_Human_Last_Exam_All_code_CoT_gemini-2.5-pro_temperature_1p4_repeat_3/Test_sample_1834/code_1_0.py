import math

# This script calculates the magnitude of the magnetic field at a specific point
# due to two infinite current-carrying wires. The final answer is symbolic.

# Define the point of interest P(x, y, z)
x, y, z = 1, -1, 0

print(f"Calculating the magnetic field at point P({x}, {y}, {z}).\n")

print("--- Step 1: Analyze the magnetic field from Wire 1 (along x-axis) ---")
# The formula for the magnetic field B from an infinite wire is B = (μ₀ * I) / (2 * π * r)
# Wire 1 has current I in the +x direction.

# The perpendicular distance r1 from the x-axis to P is its distance in the yz-plane.
r1 = math.sqrt(y**2 + z**2)
print(f"The perpendicular distance from Wire 1 to P is r1 = sqrt({y}^2 + {z}^2) = {r1}")

# By the right-hand rule (thumb in +x direction), the magnetic field at P(1, -1, 0)
# (which is below the wire) points in the negative z-direction.
print("Direction of the magnetic field B1 is in the -z direction.")
print(f"The vector B1 is therefore: - (μ₀ * I / (2 * π * {r1})) k̂")

print("\n--- Step 2: Analyze the magnetic field from Wire 2 (along y-axis) ---")
# Wire 2 has current I in the +y direction.

# The perpendicular distance r2 from the y-axis to P is its distance in the xz-plane.
r2 = math.sqrt(x**2 + z**2)
print(f"The perpendicular distance from Wire 2 to P is r2 = sqrt({x}^2 + {z}^2) = {r2}")

# By the right-hand rule (thumb in +y direction), the magnetic field at P(1, -1, 0)
# (which is to the right of the wire) points in the negative z-direction.
print("Direction of the magnetic field B2 is in the -z direction.")
print(f"The vector B2 is therefore: - (μ₀ * I / (2 * π * {r2})) k̂")

print("\n--- Step 3: Sum the vectors to find the total magnetic field ---")
# The total magnetic field is B_total = B1 + B2.
print("B_total = B1 + B2")
print(f"B_total = [- (μ₀ * I / (2 * π * {r1})) k̂] + [- (μ₀ * I / (2 * π * {r2})) k̂]")
print("Since both vectors point in the same direction, we add their magnitudes:")
print(f"B_total = - (μ₀ * I / (2 * π * {r1}) + μ₀ * I / (2 * π * {r2})) k̂")
print(f"B_total = - (μ₀ * I / (2 * π) + μ₀ * I / (2 * π)) k̂")
print("B_total = - (2 * μ₀ * I / (2 * π)) k̂")
print("B_total = - (μ₀ * I / π) k̂")

print("\n--- Step 4: Calculate the magnitude of the total magnetic field ---")
# The magnitude of B_total is the absolute value of its component.
print("|B_total| = |- (μ₀ * I / π)|")
print("The final equation for the magnitude is:")
print("Magnitude = μ₀ * I / π")
