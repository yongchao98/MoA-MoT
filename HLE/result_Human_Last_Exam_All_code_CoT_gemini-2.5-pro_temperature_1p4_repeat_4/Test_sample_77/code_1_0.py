import math

# --- Simulation Parameters ---
# Let's define the positions of the objects in our scene.
# Assume the inflow emitter sphere is centered at a height of 8.0 meters.
inflow_height = 8.0
# Assume the obstacle plane is at a height of 1.5 meters.
obstacle_height = 1.5
# Acceleration due to gravity is a constant.
gravity = 9.8

# --- Physics Calculation ---
# First, we calculate the vertical distance the fluid needs to fall.
fall_distance = inflow_height - obstacle_height

# Next, we use the kinematic equation for free fall to find the time it takes.
# The equation is: distance = 0.5 * gravity * time^2
# We can rearrange it to solve for time: time = sqrt(2 * distance / gravity)
time_to_impact = math.sqrt((2 * fall_distance) / gravity)

# --- Output Results ---
print("--- Fluid Simulation: Time to Initial Impact ---")
print(f"Inflow object height: {inflow_height} m")
print(f"Obstacle plane height: {obstacle_height} m")
print(f"Fall distance (d): {fall_distance:.2f} m")
print(f"Gravity (g): {gravity} m/s^2")
print("-" * 45)
print("The final equation to calculate time (t) is: t = sqrt(2 * d / g)")
print("Here is the equation with each number filled in:")
# The following line outputs the numbers used in the final equation.
print(f"t = sqrt(2 * {fall_distance:.2f} / {gravity})")
print(f"\nCalculated time to impact: {time_to_impact:.2f} seconds")
