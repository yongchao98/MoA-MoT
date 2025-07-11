from fractions import Fraction

# --- Problem Setup ---
# Let L be the width of the river.
# Let v_0 be the maximum flow velocity at the center.
# Let v be the boat's speed relative to the water, perpendicular to the flow.
# We set up a coordinate system:
# - The y-axis points from the starting bank (y=0) to the opposite bank (y=L).
# - The x-axis points downstream along the starting bank. The boat starts at (0,0).

# The river's velocity v_x at a distance y from the bank is given by a triangular profile:
# v_x(y) = (2*v_0/L) * y      for 0 <= y <= L/2
# v_x(y) = (2*v_0/L) * (L-y)  for L/2 <= y <= L

# --- Goal ---
# Calculate the total downstream distance drifted by the boat.

print("This script calculates the final downstream distance of the boat.")
print("The calculation is performed symbolically to derive the final formula.\n")

# --- Step-by-step Calculation ---
# The downstream drift `dx` during a small time `dt` is `dx = v_x(y) * dt`.
# The cross-stream distance traveled `dy` is `dy = v * dt`, so `dt = dy / v`.
# By substitution, `dx = v_x(y) * (dy / v)`.
# The total drift is the integral of `dx` over the boat's path in the y-direction.

# Define symbolic variables for printing the formula
L, v0, v = "L", "v_0", "v"
base_term = f"({v0}*{L}/{v})"

print("--- 1. Calculating Outbound Drift (from y=0 to y=3L/4) ---")

# The integral for drift must be split at y=L/2 because the velocity profile changes.
# Part 1: Drift from y=0 to y=L/2.
# Integral[ (2*v_0*y)/(L*v) dy ] from 0 to L/2 = (v_0*L)/(4*v)
coeff_part1 = Fraction(1, 4)
print(f"Drift during the first part of the outbound trip (y=0 to L/2) is: {coeff_part1} * {base_term}")

# Part 2: Drift from y=L/2 to y=3L/4.
# Integral[ (2*v_0*(L-y))/(L*v) dy ] from L/2 to 3L/4 = (3*v_0*L)/(16*v)
coeff_part2 = Fraction(3, 16)
print(f"Drift during the second part of the outbound trip (y=L/2 to 3L/4) is: {coeff_part2} * {base_term}")

# Total outbound drift is the sum of the two parts.
coeff_outbound = coeff_part1 + coeff_part2
print(f"Total outbound drift = ({coeff_part1} + {coeff_part2}) * {base_term} = {coeff_outbound} * {base_term}\n")

print("--- 2. Calculating Inbound Drift (from y=3L/4 to y=0) ---")

# The boat's velocity across the river is now in the -y direction, but the time taken
# to cross any segment dy is the same. The river velocity v_x(y) is also unchanged.
# Therefore, the drift during the inbound trip is the same as the outbound trip.
coeff_inbound = coeff_outbound
print(f"Due to symmetry, the inbound drift is equal to the outbound drift.")
print(f"Total inbound drift = {coeff_inbound} * {base_term}\n")

print("--- 3. Calculating Total Distance ---")
# The total distance from the origin is the sum of the outbound and inbound drifts.
total_coeff = coeff_outbound + coeff_inbound
print(f"Total Distance = Outbound Drift + Inbound Drift")
print(f"Total Distance = ({coeff_outbound} * {base_term}) + ({coeff_inbound} * {base_term})")
print(f"Total Distance = {total_coeff} * {base_term}\n")

print("--- Final Answer ---")
print("The total distance between the boat's returning position and its original starting point is given by the final equation:")
# The following line prints the final derived equation with all its numbers and variables.
print(f"Distance = ({total_coeff.numerator} * {v0} * {L}) / ({total_coeff.denominator} * {v})")