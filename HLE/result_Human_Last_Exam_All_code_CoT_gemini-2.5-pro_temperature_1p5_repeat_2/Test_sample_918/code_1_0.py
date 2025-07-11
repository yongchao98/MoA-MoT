import math

# Define the given parameters
tunnel_length_m = 2.0
moth_speed_air_mpm = 5.675  # m/min
wind_speed_mpm = 5.0      # m/min, moving east
num_leds = 80
led_delay_s = 0.3         # seconds

# --- Step 1: Calculate the moth's effective speed relative to the ground ---
# The moth flies west, against the eastward wind.
moth_speed_ground_mpm = moth_speed_air_mpm - wind_speed_mpm
print(f"Step 1: The moth's effective speed towards the west is {moth_speed_air_mpm} m/min - {wind_speed_mpm} m/min = {moth_speed_ground_mpm} m/min.")
print("-" * 20)

# --- Step 2: Calculate the time for the moth to reach the halfway point ---
# The halfway point is half the tunnel's length.
halfway_distance_m = tunnel_length_m / 2.0
time_to_halfway_min = halfway_distance_m / moth_speed_ground_mpm
print(f"Step 2: The time for the moth to reach the halfway point ({halfway_distance_m}m) is {halfway_distance_m}m / {moth_speed_ground_mpm} m/min = {time_to_halfway_min:.4f} minutes.")
print("-" * 20)

# --- Step 3: Calculate the light signal propagation time ---
# The signal travels across 79 gaps between the 80 LEDs.
light_propagation_time_s = (num_leds - 1) * led_delay_s
# Convert this time to minutes to match other units.
light_propagation_time_min = light_propagation_time_s / 60.0
print(f"Step 3: The time for the light signal to travel from the western to the eastern end is ({num_leds} - 1) * {led_delay_s}s = {light_propagation_time_s:.1f} seconds.")
print(f"   This is equal to {light_propagation_time_s:.1f}s / 60s/min = {light_propagation_time_min:.4f} minutes.")
print("-" * 20)

# --- Step 4: Calculate the total elapsed time ---
total_time_min = time_to_halfway_min + light_propagation_time_min
print(f"Step 4: The total time until the easternmost LED blinks is the sum of the time to halfway and the light propagation time.")
print(f"   Total time = {time_to_halfway_min:.4f} min + {light_propagation_time_min:.4f} min = {total_time_min:.4f} minutes.")
print("-" * 20)

# --- Step 5: Calculate the moth's final displacement ---
# Displacement = Speed * Time
final_displacement_m = moth_speed_ground_mpm * total_time_min

# We are asked to output each number in the final equation.
# Final displacement = moth_speed_ground * (time_to_halfway + light_propagation_time)
print("Step 5: The final displacement of the moth is calculated as its ground speed multiplied by the total time.")
print("The final equation is:")
# The calculation can be simplified to 1m + displacement during light propagation
# But to show the full equation as requested:
# Displacement = 0.675 * (1/0.675 + 0.395)
print(f"Displacement = {moth_speed_ground_mpm} m/min * ({time_to_halfway_min:.4f} min + {light_propagation_time_min:.4f} min)")
print(f"Displacement = {moth_speed_ground_mpm} m/min * {total_time_min:.4f} min")
print(f"Final Displacement = {final_displacement_m:.2f} m")

<<<D>>>