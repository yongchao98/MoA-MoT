import math

# --- Symbolic Representation ---
# We will use strings to represent the symbolic constants.
mu_0 = "μ₀"
I = "I"
pi = "π"

# --- Given Point ---
x, y, z = 1, -1, 0

print(f"Calculating the magnetic field at point P({x}, {y}, {z}).\n")

# --- Step 1: Contribution from Wire 1 (on x-axis) ---
print("Step 1: Analyze the magnetic field from Wire 1 (on x-axis).")
# The perpendicular distance from a point (x, y, z) to the x-axis is sqrt(y^2 + z^2).
r1 = math.sqrt(y**2 + z**2)
print(f"The perpendicular distance (r1) from Wire 1 to P is sqrt({y}^2 + {z}^2) = {int(r1)}.")

# The magnitude of the magnetic field from Wire 1 is B1 = (μ₀ * I) / (2 * π * r1).
print(f"The magnitude of the field B1 is given by the formula (μ₀ * I) / (2 * π * r1).")
print(f"Substituting r1 = {int(r1)}, we get: B1 = ({mu_0} * {I}) / (2 * {pi} * {int(r1)}).")

# Determine the direction of B1 using the right-hand rule.
# Current is in +x. The point is at y = -1. Fingers curl into the page (-z direction).
print("By the right-hand rule, the direction of B1 at point P is in the -z direction.")
print("-" * 40)

# --- Step 2: Contribution from Wire 2 (on y-axis) ---
print("Step 2: Analyze the magnetic field from Wire 2 (on y-axis).")
# The perpendicular distance from a point (x, y, z) to the y-axis is sqrt(x^2 + z^2).
r2 = math.sqrt(x**2 + z**2)
print(f"The perpendicular distance (r2) from Wire 2 to P is sqrt({x}^2 + {z}^2) = {int(r2)}.")

# The magnitude of the magnetic field from Wire 2 is B2 = (μ₀ * I) / (2 * π * r2).
print(f"The magnitude of the field B2 is given by the formula (μ₀ * I) / (2 * π * r2).")
print(f"Substituting r2 = {int(r2)}, we get: B2 = ({mu_0} * {I}) / (2 * {pi} * {int(r2)}).")

# Determine the direction of B2 using the right-hand rule.
# Current is in +y. The point is at x = 1. Fingers curl into the page (-z direction).
print("By the right-hand rule, the direction of B2 at point P is also in the -z direction.")
print("-" * 40)

# --- Step 3: Total Magnetic Field ---
print("Step 3: Calculate the total magnetic field magnitude.")
print("Since both field vectors B1 and B2 point in the same direction (-z), we can add their magnitudes to find the total magnitude.")
print("Magnitude |B_total| = B1 + B2")
print(f"|B_total| = [({mu_0} * {I}) / (2 * {pi} * {int(r1)})] + [({mu_0} * {I}) / (2 * {pi} * {int(r2)})]")
print(f"Since r1 = r2 = {int(r1)}, this simplifies to:")
print(f"|B_total| = 2 * (({mu_0} * {I}) / (2 * {pi} * {int(r1)}))")
print(f"|B_total| = ({mu_0} * {I}) / ({pi} * {int(r1)})")

# Final result
final_magnitude_expression = f"({mu_0} * {I}) / {pi}"
print("\nFinal Answer:")
print(f"The magnitude of the magnetic field at ({x}, {y}, {z}) is {final_magnitude_expression}.")
