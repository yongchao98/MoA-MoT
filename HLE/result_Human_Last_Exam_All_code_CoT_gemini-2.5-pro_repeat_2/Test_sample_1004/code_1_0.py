import sys

# This script calculates the downstream drift of a boat in a river with a parabolic flow profile.

# --- Problem Parameters (user can modify these) ---
L = 100.0   # Width of the river in meters
v_0 = 5.0   # Maximum flow velocity at the center in m/s
v = 4.0     # Boat's speed perpendicular to the flow (relative to water) in m/s

# --- Derivation Summary ---
# The river flow velocity is u(y) = (4*v_0/L^2) * y * (L-y).
# The boat travels out from y=0 to y=3L/4, then returns from y=3L/4 to y=0.
# The total downstream drift is the sum of the drift from both legs of the journey.
# Drift = integral of u(y) dt. Since dt can be written as dy/v (in magnitude), the
# drift calculation for both legs turns out to be identical.
# Total Drift = 2 * integral from 0 to 3L/4 of [u(y)/v] dy
# After solving the integral, the final formula is:
# Total Drift = (9 * v_0 * L) / (8 * v)

# --- Calculation ---
# The constant factors (numbers) in the final equation are 9 and 8.
numerator_factor = 9
denominator_factor = 8

# Check for division by zero to prevent runtime errors.
if v == 0:
    print("Error: The boat's relative speed 'v' cannot be zero.", file=sys.stderr)
    sys.exit(1)

# Apply the final derived formula to calculate the total drift distance.
total_drift = (numerator_factor * v_0 * L) / (denominator_factor * v)

# --- Output ---
print("This script calculates the total downstream drift for the described boat journey.")
print("\n--- Input Parameters ---")
print(f"River width (L): {L} m")
print(f"Max flow speed (v_0): {v_0} m/s")
print(f"Boat's relative speed (v): {v} m/s")

print("\n--- Result ---")
# As requested, printing the numbers in the final equation explicitly.
print(f"The calculation is based on the formula: ( {numerator_factor} * v_0 * L ) / ( {denominator_factor} * v )")
print(f"The total distance from the start point is: {total_drift:.2f} meters")