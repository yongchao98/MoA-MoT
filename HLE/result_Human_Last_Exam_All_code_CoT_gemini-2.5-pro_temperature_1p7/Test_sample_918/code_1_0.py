import math

# Define the given parameters
tunnel_length_m = 2.0
wind_speed_mpm = 5.0
moth_speed_in_air_mpm = 5.675
num_leds = 80
led_delay_s = 0.3
seconds_per_minute = 60.0

# Step 1: Calculate the moth's effective speed relative to the ground.
# The moth flies west (against the wind), and the wind moves east.
moth_ground_speed_mpm = moth_speed_in_air_mpm - wind_speed_mpm

# Step 2: Calculate the time it takes for the moth to reach the halfway point.
# The moth starts at the eastern end and needs to travel to the center.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm

# Step 3: Calculate the time for the LED signal to travel from west to east.
# The signal travels across 80 LEDs, which means there are 79 intervals.
num_intervals = num_leds - 1
signal_travel_time_s = num_intervals * led_delay_s
signal_travel_time_min = signal_travel_time_s / seconds_per_minute

# Step 4: Calculate the total time elapsed from the moth's start.
total_time_min = time_to_halfway_min + signal_travel_time_min

# Step 5: Calculate the moth's total displacement during this total time.
total_displacement_m = moth_ground_speed_mpm * total_time_min

# Print the step-by-step calculation
print("Step 1: Calculate the moth's effective ground speed.")
print(f"Moth's speed - Wind speed = {moth_speed_in_air_mpm} m/min - {wind_speed_mpm} m/min = {moth_ground_speed_mpm:.3f} m/min")
print("-" * 30)

print("Step 2: Calculate the time to reach the halfway point.")
print(f"Time to halfway = Distance / Speed = {distance_to_halfway_m:.1f} m / {moth_ground_speed_mpm:.3f} m/min = {time_to_halfway_min:.4f} min")
print("-" * 30)

print("Step 3: Calculate the time for the LED signal to cross the tunnel.")
print(f"Signal travel time = ({num_leds} - 1) * {led_delay_s} s = {signal_travel_time_s:.1f} s")
print(f"In minutes: {signal_travel_time_s:.1f} s / {seconds_per_minute} s/min = {signal_travel_time_min:.4f} min")
print("-" * 30)

print("Step 4: Calculate the final displacement.")
print(f"Displacement = Ground Speed * (Time to Halfway + Signal Travel Time)")
print(f"Displacement = {moth_ground_speed_mpm:.3f} m/min * ({time_to_halfway_min:.4f} min + {signal_travel_time_min:.4f} min)")
print(f"Displacement = {moth_ground_speed_mpm:.3f} m/min * {total_time_min:.4f} min")
print(f"Final Displacement = {total_displacement_m:.2f} m")

<<<D>>>