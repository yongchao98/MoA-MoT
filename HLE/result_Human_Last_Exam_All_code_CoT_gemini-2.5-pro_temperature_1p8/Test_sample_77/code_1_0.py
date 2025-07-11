import math

# Step 1: Define scene parameters and constants.
# Let's assume the inflow object is 4 meters above the obstacle plane.
distance = 4.0  # meters
# Acceleration due to gravity.
gravity = 9.8   # meters/second^2
# The initial vertical velocity is 0 as the fluid just starts to fall.
initial_velocity = 0.0

# Step 2 & 3: Calculate the time to impact using the simplified kinematic formula.
# The formula is t = sqrt(2 * d / a)
time_to_impact = math.sqrt(2 * distance / gravity)

# Step 4: Print the results, including the numbers in the final equation.
print(f"Calculating the time for a fluid particle to fall {distance}m under {gravity} m/s^2 gravity.")
print("The governing equation for time (t) is: t = sqrt(2 * distance / gravity)")
print(f"Plugging in the numbers: t = sqrt(2 * {distance} / {gravity})")
print(f"Time to impact the obstacle: {time_to_impact:.4f} seconds")

# The final answer format is not directly related to the textual explanation,
# but it will be the numerical result of the calculation.