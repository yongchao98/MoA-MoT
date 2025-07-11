import math

# Define the given constants
moth_airspeed = 5.675  # m/min
wind_speed = 5.0       # m/min
tunnel_length = 2.0    # m
num_lights = 80
delay_between_lights_s = 0.3 # seconds

# Step 1: Calculate the moth's speed relative to the ground.
# The moth flies against the wind.
moth_ground_speed = moth_airspeed - wind_speed

# Step 2: Calculate the time for the moth to reach the halfway point.
# The halfway point is 1m from the starting eastern end.
distance_to_halfway = tunnel_length / 2.0
time_to_halfway_min = distance_to_halfway / moth_ground_speed

# Step 3: Calculate the time for the LED signal to travel across the tunnel.
# There are (num_lights - 1) intervals between the lights.
signal_travel_time_s = (num_lights - 1) * delay_between_lights_s
# Convert seconds to minutes for unit consistency.
signal_travel_time_min = signal_travel_time_s / 60.0

# Step 4: Calculate the total time until the easternmost LED blinks.
# This is the time to reach the halfway point plus the signal travel time.
total_time_min = time_to_halfway_min + signal_travel_time_min

# Step 5: Calculate the moth's total displacement from its starting point.
# Displacement = speed * total time
displacement_m = moth_ground_speed * total_time_min

# Print the final equation with all the numbers.
# The displacement can be expanded as:
# displacement = moth_ground_speed * (time_to_halfway + signal_travel_time_min)
# displacement = moth_ground_speed * (distance_to_halfway / moth_ground_speed + signal_travel_time_min)
# displacement = distance_to_halfway + moth_ground_speed * signal_travel_time_min
# This shows that the moth travels 1m, and then travels for an additional duration equal to the signal propagation time.

final_calc_part1 = distance_to_halfway
final_calc_part2 = moth_ground_speed * signal_travel_time_min

print(f"The moth's displacement is the distance to the halfway point plus the distance traveled during the LED signal propagation.")
print(f"Final Displacement = {final_calc_part1:.2f}m + ({moth_airspeed}m/min - {wind_speed}m/min) * (({num_lights} - 1) * {delay_between_lights_s}s / 60s/min)")
print(f"Final Displacement = {final_calc_part1:.2f}m + {moth_ground_speed:.3f}m/min * {signal_travel_time_min:.3f}min")
print(f"Final Displacement = {final_calc_part1:.2f}m + {final_calc_part2:.4f}m = {displacement_m:.4f}m")
print(f"The moth's displacement is approximately {displacement_m:.2f}m.")
<<<D>>>