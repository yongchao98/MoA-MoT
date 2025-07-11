import math

# Define the point of interest P(x,y,z)
x, y, z = 1, -1, 0

# --- Calculation for Wire 1 (on x-axis) ---
# The perpendicular distance from P to the x-axis is sqrt(y^2 + z^2)
r1 = math.sqrt(y**2 + z**2)

# The magnitude of the magnetic field B1 is given by B1 = (mu_0 * I) / (2 * pi * r1)
# We will represent this as a coefficient multiplied by (mu_0 * I)
# The coefficient is 1 / (2 * pi * r1)
coeff1 = 1 / (2 * math.pi * r1)

# --- Calculation for Wire 2 (on y-axis) ---
# The perpendicular distance from P to the y-axis is sqrt(x^2 + z^2)
r2 = math.sqrt(x**2 + z**2)

# The magnitude of the magnetic field B2 is given by B2 = (mu_0 * I) / (2 * pi * r2)
# The coefficient is 1 / (2 * pi * r2)
coeff2 = 1 / (2 * math.pi * r2)

# --- Total Magnetic Field Magnitude ---
# As determined by the right-hand rule, both B1 and B2 point in the -z direction.
# Therefore, the magnitude of the total magnetic field is the sum of the individual magnitudes.
total_magnitude_coeff = coeff1 + coeff2

print("Calculation Steps:")
print(f"Point of interest P = ({x}, {y}, {z})")
print("\n--- For Wire 1 (on x-axis) ---")
print(f"Distance from P to Wire 1, r1 = {r1}")
print(f"Magnitude B1 = (μ₀ * I) / (2 * π * {r1})")

print("\n--- For Wire 2 (on y-axis) ---")
print(f"Distance from P to Wire 2, r2 = {r2}")
print(f"Magnitude B2 = (μ₀ * I) / (2 * π * {r2})")

print("\n--- Total Magnitude ---")
print("The magnetic fields from both wires point in the same direction (-z).")
print("So, the total magnitude is the sum of the individual magnitudes:")
print(f"|B_total| = B1 + B2")
print(f"|B_total| = (μ₀ * I) / (2 * π * {r1}) + (μ₀ * I) / (2 * π * {r2})")
print(f"|B_total| = 2 * (μ₀ * I) / (2 * π)")
print(f"|B_total| = μ₀ * I / π")