import math

# Define the initial parameters from the problem description
tunnel_length_m = 2.0
wind_speed_mpmin = 5.0
moth_airspeed_mpmin = 5.675
num_leds = 80
led_blink_delay_s = 0.3

# --- Step 1 & 2: Calculate Moth's Ground Speed in m/s ---
# Convert speeds from meters per minute to meters per second
wind_speed_mps = wind_speed_mpmin / 60.0
moth_airspeed_mps = moth_airspeed_mpmin / 60.0

# The moth flies against the wind, so we subtract the wind speed from the moth's airspeed.
moth_ground_speed_mps = moth_airspeed_mps - wind_speed_mps

# --- Step 3: Calculate Time to Reach Halfway Point ---
# The blinking sequence starts when the moth is halfway through the tunnel.
halfway_distance_m = tunnel_length_m / 2.0

# Time = Distance / Speed
time_to_halfway_s = halfway_distance_m / moth_ground_speed_mps

# --- Step 4: Calculate LED Signal Propagation Time ---
# The signal travels across the (num_leds - 1) gaps between the lights.
num_gaps = num_leds - 1
led_propagation_time_s = num_gaps * led_blink_delay_s

# --- Step 5: Calculate Total Time ---
# The total time is the time for the moth to reach the halfway point plus the time for the blink signal to travel.
total_time_s = time_to_halfway_s + led_propagation_time_s

# --- Step 6: Calculate Final Displacement ---
# The final displacement is the moth's ground speed multiplied by the total time.
final_displacement_m = moth_ground_speed_mps * total_time_s

# --- Output the results ---
print("This script calculates the moth's displacement based on the provided conditions.")
print("-" * 50)
print(f"1. Moth's speed relative to the ground: {moth_airspeed_mpmin} m/min - {wind_speed_mpmin} m/min = {moth_ground_speed_mps * 60:.3f} m/min, or {moth_ground_speed_mps:.5f} m/s.")
print(f"2. Time for moth to travel to halfway point (1.0m): {halfway_distance_m:.1f} m / {moth_ground_speed_mps:.5f} m/s = {time_to_halfway_s:.2f} s.")
print(f"3. Time for LED blink to travel from west to east end: ({num_leds} - 1) * {led_blink_delay_s} s = {led_propagation_time_s:.2f} s.")
print(f"4. Total time elapsed: {time_to_halfway_s:.2f} s + {led_propagation_time_s:.2f} s = {total_time_s:.2f} s.")
print("-" * 50)
print("Final displacement calculation:")
print(f"Displacement = Moth's ground speed * Total time")
print(f"Displacement = {moth_ground_speed_mps:.5f} m/s * {total_time_s:.2f} s")
print(f"Final Displacement: {final_displacement_m:.2f} m")

<<<D>>>