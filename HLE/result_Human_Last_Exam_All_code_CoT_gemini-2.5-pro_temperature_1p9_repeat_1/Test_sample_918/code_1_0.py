import math

# Define the given variables
tunnel_length_m = 2.0
moth_airspeed_mpm = 5.675  # meters per minute
wind_speed_mpm = 5.0      # meters per minute
num_leds = 80
led_delay_s = 0.3         # seconds
led_blink_interval_s = 2.0

# Step 1: Calculate the moth's effective speed relative to the ground.
# The moth flies west, against the eastbound wind.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
print(f"The moth's effective speed against the wind is {moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min = {moth_ground_speed_mpm} m/min.")
print("-" * 20)

# Step 2: Calculate the time it takes for the moth to reach the halfway point.
# The moth starts at the eastern end (2m mark) and the halfway point is at the 1m mark.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm
print(f"The LED sequence starts when the moth is halfway, so it must first travel {distance_to_halfway_m} m.")
print(f"Time to reach halfway point = {distance_to_halfway_m} m / {moth_ground_speed_mpm} m/min = {time_to_halfway_min:.4f} min.")
print("-" * 20)

# Step 3: Calculate the time it takes for the LED signal to propagate to the end.
# The signal travels across (num_leds - 1) gaps.
propagation_time_s = (num_leds - 1) * led_delay_s
propagation_time_min = propagation_time_s / 60.0
print(f"The blinking signal starts at the west end and must travel to the east end.")
print(f"Total propagation time = ({num_leds} - 1) * {led_delay_s} s = {propagation_time_s} s.")
print(f"In minutes, this is {propagation_time_s} s / 60 s/min = {propagation_time_min} min.")
print("-" * 20)

# Step 4: Calculate the total time the moth is in flight.
total_time_min = time_to_halfway_min + propagation_time_min
print(f"The moth's total flight time is the sum of the time to reach the halfway point and the signal propagation time.")
print(f"Total time = {time_to_halfway_min:.4f} min + {propagation_time_min} min = {total_time_min:.4f} min.")
print("-" * 20)

# Step 5: Calculate the moth's total displacement.
# Displacement = speed * time
displacement_m = moth_ground_speed_mpm * total_time_min

print("The moth's total displacement is calculated by multiplying its ground speed by the total flight time.")
# Final equation output as requested
print(f"Displacement = {moth_ground_speed_mpm} m/min * ({time_to_halfway_min:.4f} min + {propagation_time_min} min)")
print(f"Displacement = {moth_ground_speed_mpm} m/min * {total_time_min:.4f} min")
print(f"Final Displacement = {displacement_m:.2f} m")

<<<D>>>