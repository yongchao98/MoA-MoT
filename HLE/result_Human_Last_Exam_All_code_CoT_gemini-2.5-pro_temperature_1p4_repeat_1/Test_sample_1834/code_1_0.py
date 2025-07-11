import math

# This script calculates the magnitude of the magnetic field at a specific point
# due to two perpendicular infinite current-carrying wires.

# --- Define the point of interest (x, y, z) ---
x, y, z = 1, -1, 0

# --- Step 1: Analyze Wire 1 (on x-axis) ---
# The perpendicular distance, r1, from the point to the x-axis is sqrt(y^2 + z^2).
r1 = math.sqrt(y**2 + z**2)
# The magnitude of the field from wire 1 is B1 = (mu_0 * I) / (2 * pi * r1).

# --- Step 2: Analyze Wire 2 (on y-axis) ---
# The perpendicular distance, r2, from the point to the y-axis is sqrt(x^2 + z^2).
r2 = math.sqrt(x**2 + z**2)
# The magnitude of the field from wire 2 is B2 = (mu_0 * I) / (2 * pi * r2).

# --- Step 3: Calculate Total Magnetic Field Magnitude ---
# At the given point, the magnetic fields from both wires point in the same direction (-z).
# Therefore, the magnitude of the total field is the sum of their individual magnitudes.
# B_total = B1 + B2

# We will now construct and print the final equation, including the calculated numbers.
print("The formula for the total magnetic field magnitude is B_total = B1 + B2.")
print("B_total = (mu_0 * I) / (2 * pi * r1) + (mu_0 * I) / (2 * pi * r2)")
print("\nFirst, we calculate the distances r1 and r2:")
print(f"r1 = sqrt(({y})^2 + ({z})^2) = {r1}")
print(f"r2 = sqrt(({x})^2 + ({z})^2) = {r2}")
print("\nSubstituting these numbers into the equation gives:")
# We use .1f to format the numbers as floats for clarity.
print(f"B_total = (mu_0 * I) / (2 * pi * {r1:.1f}) + (mu_0 * I) / (2 * pi * {r2:.1f})")

# Simplifying the expression:
# B_total = (mu_0 * I / (2 * pi)) * (1/r1 + 1/r2)
# Since r1=1 and r2=1, this is (mu_0 * I / (2 * pi)) * (1 + 1) = (mu_0 * I / (2 * pi)) * 2
# which simplifies to (mu_0 * I) / pi.
# So, the final coefficient multiplying (mu_0 * I) / pi is 1.0.
final_coefficient = (1.0 / (2 * math.pi * r1) + 1.0 / (2 * math.pi * r2)) * math.pi

print("\nSimplifying this expression leads to the final answer:")
print(f"B_total = {final_coefficient:.1f} * (mu_0 * I) / pi")