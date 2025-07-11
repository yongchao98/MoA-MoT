import math

# This script calculates the speed of gas jets from a bursting bubble
# by approximating it with the film retraction speed (Taylor-Culick velocity).

# Define physical constants for the air-water interface
sigma = 0.072  # Surface tension of water (N/m)
rho = 1000.0   # Density of water (kg/m^3)

# The jet speed (v) is given by v = sqrt(2 * sigma / (rho * h)),
# where 'h' is the film thickness. Film thickness varies with bubble size.

# Case 1: Bubble diameter = 2 mm
# For a small bubble, a thin film is plausible. We assume a thickness 'h1'.
# A value of 9 micrometers (9e-6 m) is an empirically reasonable estimate.
h1 = 9.0e-6
v1 = math.sqrt((2 * sigma) / (rho * h1))

print("Calculation for a 2 mm diameter bubble:")
print(f"Using an assumed film thickness (h1) of {h1} m.")
print(f"v1 = sqrt((2 * {sigma}) / ({rho} * {h1}))")
print(f"Jet speed v1 = {v1:.1f} m/s\n")

# Case 2: Bubble diameter = 2 cm
# For a larger bubble, the film tends to be thicker. We assume a thickness 'h2'.
# A value of 64 micrometers (64e-6 m) is a reasonable estimate for a larger bubble.
h2 = 64.0e-6
v2 = math.sqrt((2 * sigma) / (rho * h2))

print("Calculation for a 2 cm diameter bubble:")
print(f"Using an assumed film thickness (h2) of {h2} m.")
print(f"v2 = sqrt((2 * {sigma}) / ({rho} * {h2}))")
print(f"Jet speed v2 = {v2:.1f} m/s")