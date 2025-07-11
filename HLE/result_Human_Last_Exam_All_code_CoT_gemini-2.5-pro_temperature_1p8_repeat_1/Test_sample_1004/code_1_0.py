from fractions import Fraction

# --- Problem Setup ---
# Let the river width be L.
# Let the starting bank be at y = 0 and the opposite bank at y = L.
# The boat starts at the origin (x=0, y=0).
# The boat's speed relative to the water is v, always perpendicular to the current.
# The river's maximum flow velocity is v0 at the center (y = L/2).

# Symbolic variables for printing the final equation
L_str = "L"
v0_str = "v0"
v_str = "v"

# --- Step 1: Define River Velocity v_river(y) ---
# For 0 <= y <= L/2, v_river(y) = (2 * v0 / L) * y
# For L/2 < y <= L, v_river(y) = (2 * v0 / L) * (L - y)
# The drift dx at a given y is dx = v_river(y) * dt.
# The time to cross dy is dt = dy / v.
# So, the total drift for a one-way trip is the integral of (v_river(y) / v) dy.

print("This script calculates the total downstream drift of the boat.")
print("The calculation involves integrating the river's velocity profile.\n")

# --- Step 2: Calculate Drift for the Outbound Trip (y=0 to y=3L/4) ---
# The total drift is (1/v) * Integral(v_river(y) dy) from y=0 to y=3L/4.
# We calculate the coefficient of the (v0 * L) term from the integral.

# Part 1 of the integral: from y=0 to y=L/2
# Integral of (2/L)*y dy = y^2/L. Evaluated from 0 to L/2: (L/2)^2/L = L/4.
# The coefficient of (v0*L) is 1/4.
integral_coeff_1 = Fraction(1, 4)

# Part 2 of the integral: from y=L/2 to y=3L/4
# Integral of (2/L)*(L-y) dy = (2/L)*(L*y - y^2/2). Evaluated from L/2 to 3L/4.
# Value at 3L/4: (2/L)*(L*(3L/4) - (3L/4)^2/2) = (2/L)*(3L^2/4 - 9L^2/32) = (2/L)*(15L^2/32) = 15L/16.
# Value at L/2: (2/L)*(L*(L/2) - (L/2)^2/2) = (2/L)*(L^2/2 - L^2/8) = (2/L)*(3L^2/8) = 3L/4 = 12L/16.
# Result of integral part 2 is (15L/16 - 12L/16) = 3L/16.
# The coefficient of (v0*L) is 3/16.
integral_coeff_2 = Fraction(3, 16)

# The total coefficient for the one-way trip integral
oneway_integral_coeff = integral_coeff_1 + integral_coeff_2

# The drift for one trip is (integral_coeff) * v0 * L / v.
oneway_drift_coeff = oneway_integral_coeff

print(f"The integral of the velocity profile for a one-way trip from y=0 to y=3L/4 gives a term proportional to {oneway_drift_coeff.numerator}/{oneway_drift_coeff.denominator} * {v0_str} * {L_str}.")

# --- Step 3: Calculate Total Drift ---
# The boat travels from y=0 to y=3L/4 and then back to y=0.
# The time spent at each y-coordinate is the same for both legs of the journey.
# Therefore, the drift on the return trip is the same as the drift on the outbound trip.
# Total Drift = 2 * (Outbound Drift)
total_drift_coeff = 2 * oneway_drift_coeff

# --- Step 4: Display the Final Answer ---
# The problem asks to output each number in the final equation.
numerator = total_drift_coeff.numerator
denominator = total_drift_coeff.denominator

print("\nThe total drift is twice the one-way drift.")
print("\nFinal Result:")
print(f"The total distance between the start and return points is given by the equation:")
print(f"Total Drift = ({numerator} / {denominator}) * ({v0_str} * {L_str} / {v_str})")
print(f"Total Drift = ({numerator} * {v0_str} * {L_str}) / ({denominator} * {v_str})")
