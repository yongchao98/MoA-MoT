import math

# --- Define problem constants ---
tunnel_length_m = 2.0
wind_speed_mpm = 5.0
moth_airspeed_mpm = 5.675
num_leds = 80
led_delay_s = 0.3
seconds_per_minute = 60.0

# --- Step 1: Calculate the moth's net speed relative to the ground ---
# The moth flies west, against the eastbound wind.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

# --- Step 2 & 3: Calculate the two time components ---
# Time for moth to travel to the halfway point (trigger point for lights)
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm

# Time for the light wave to propagate from the first to the last LED
time_for_wave_s = (num_leds - 1) * led_delay_s
time_for_wave_min = time_for_wave_s / seconds_per_minute

# --- Step 4: Calculate the total displacement ---
# Total displacement is the distance to halfway + distance traveled during the wave
# which is equivalent to ground speed * total time.
total_time_min = time_to_halfway_min + time_for_wave_min
total_displacement_m = moth_ground_speed_mpm * total_time_min

# --- Step 5: Print the final equation and result ---
print("The moth's total displacement is its ground speed multiplied by the total time elapsed.")
print("Total time is the time to reach the halfway point plus the time for the light wave to travel.")
print("\nFinal Equation:")
print("Displacement = (Moth Airspeed - Wind Speed) * [ (Tunnel Length / 2) / (Moth Airspeed - Wind Speed) + ((Num LEDs - 1) * LED Delay) / Seconds per Minute ]")
print(f"Displacement = ({moth_airspeed_mpm} - {wind_speed_mpm}) * [ ({tunnel_length_m} / 2) / ({moth_airspeed_mpm} - {wind_speed_mpm}) + (({num_leds} - 1) * {led_delay_s}) / {seconds_per_minute} ]")

# A simpler way to view the calculation is the sum of two distances:
dist_during_wave = moth_ground_speed_mpm * time_for_wave_min
print("\nSimplified Calculation:")
print("Displacement = Distance to halfway + Distance traveled during light wave")
print(f"Displacement = {distance_to_halfway_m} m + {moth_ground_speed_mpm} m/min * {time_for_wave_min:.3f} min")
print(f"Displacement = {distance_to_halfway_m} m + {dist_during_wave:.4f} m")

# Final result rounded to two decimal places for the answer choice.
print(f"\nThe moth's total displacement is: {total_displacement_m:.2f} m")
