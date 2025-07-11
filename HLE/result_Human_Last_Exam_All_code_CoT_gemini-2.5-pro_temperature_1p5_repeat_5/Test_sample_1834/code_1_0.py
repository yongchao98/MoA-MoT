import numpy as np

# Define the point of interest P
px, py, pz = 1, -1, 0
P = (px, py, pz)

print("This script calculates the magnitude of the magnetic field at a point due to two infinite wires.\n")

# --- Contribution from Wire 1 (on x-axis) ---
print("Step 1: Analyze the magnetic field from Wire 1 (along x-axis, current in +x direction).")
# The perpendicular distance r1 from point P(x,y,z) to the x-axis is sqrt(y^2 + z^2).
r1 = np.sqrt(py**2 + pz**2)
print(f"The perpendicular distance from point P{P} to Wire 1 is r1 = sqrt(({py})^2 + ({pz})^2) = {r1}.")

# The magnitude of the magnetic field from Wire 1 is B1 = (μ₀ * I) / (2 * π * r1).
# The direction is found using the right-hand rule. With the thumb in the +x direction,
# the fingers curl. At the point P(1, -1, 0), which is below the x-axis, the field points in the +z direction.
print("Using the right-hand rule, the field B1 at P points in the positive z-direction.")
print("The vector B1 is B1 = ( (μ₀ * I) / (2 * π * r1) ) * k̂\n")


# --- Contribution from Wire 2 (on y-axis) ---
print("Step 2: Analyze the magnetic field from Wire 2 (along y-axis, current in +y direction).")
# The perpendicular distance r2 from point P(x,y,z) to the y-axis is sqrt(x^2 + z^2).
r2 = np.sqrt(px**2 + pz**2)
print(f"The perpendicular distance from point P{P} to Wire 2 is r2 = sqrt(({px})^2 + ({pz})^2) = {r2}.")

# The magnitude of the magnetic field from Wire 2 is B2 = (μ₀ * I) / (2 * π * r2).
# The direction is found using the right-hand rule. With the thumb in the +y direction,
# the fingers curl. At the point P(1, -1, 0), which is to the right of the y-axis, the field points in the -z direction.
print("Using the right-hand rule, the field B2 at P points in the negative z-direction.")
print("The vector B2 is B2 = - ( (μ₀ * I) / (2 * π * r2) ) * k̂\n")


# --- Total Magnetic Field ---
print("Step 3: Calculate the total magnetic field by vector addition (Superposition Principle).")
print("B_total = B1 + B2")
print("B_total = [ (μ₀ * I) / (2 * π * r1) ] * k̂  -  [ (μ₀ * I) / (2 * π * r2) ] * k̂")

# The problem requests to show the numbers in the final equation.
print("\nSubstituting the calculated distances r1 and r2:")
print(f"B_total = [ (μ₀ * I) / (2 * π * {r1}) - (μ₀ * I) / (2 * π * {r2}) ] * k̂")

# Since r1 and r2 are both 1, the two terms are equal in magnitude and opposite in direction.
print("\nSince r1 = r2 = 1, the two terms cancel each other out:")
magnitude_term = "(μ₀ * I) / (2 * π)"
print(f"B_total = ( {magnitude_term} - {magnitude_term} ) * k̂ = 0")

final_magnitude = 0.0
print(f"\nThe magnitude of the total magnetic field at point {P} is {final_magnitude}.")
