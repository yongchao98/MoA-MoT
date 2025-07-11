import math

# Define the constants from the problem
tunnel_length_m = 2.0
wind_speed_m_per_min = 5.0
moth_air_speed_m_per_min = 5.675
num_leds = 80
led_delay_s = 0.3

# --- Step 1: Calculate the moth's speed relative to the ground ---
# The moth flies west, and the wind blows east, so their speeds subtract.
moth_ground_speed_m_per_min = moth_air_speed_m_per_min - wind_speed_m_per_min
# Convert speed to meters per second for consistency with other units
moth_ground_speed_m_per_s = moth_ground_speed_m_per_min / 60.0

print(f"The moth's air speed is {moth_air_speed_m_per_min} m/min.")
print(f"The wind speed is {wind_speed_m_per_min} m/min.")
print("Calculating the moth's actual speed relative to the ground (ground speed):")
print(f"Ground Speed = {moth_air_speed_m_per_min} m/min - {wind_speed_m_per_min} m/min = {moth_ground_speed_m_per_min} m/min\n")

# --- Step 2: Determine the time it takes for the LED signal to propagate ---
# The signal starts at the 1st LED and ends at the 80th, crossing 79 gaps.
num_gaps = num_leds - 1
led_propagation_time_s = num_gaps * led_delay_s

print("The LED blinking sequence starts when the moth is at the halfway point (1m from the start).")
print("We need to find how much farther the moth travels while the LED signal propagates from the west end to the east end.")
print("Calculating the LED signal propagation time:")
print(f"Propagation Time = (Number of LEDs - 1) * Delay per LED")
print(f"Propagation Time = ({num_leds} - 1) * {led_delay_s} s = {led_propagation_time_s:.1f} s\n")

# --- Step 3: Calculate the additional distance the moth travels ---
# Additional distance = ground speed in m/s * LED propagation time in seconds
additional_distance_m = moth_ground_speed_m_per_s * led_propagation_time_s

print("Calculating the additional distance the moth flies during this time:")
print(f"Additional Distance = Ground Speed * Propagation Time")
# Show the calculation with numbers
print(f"Additional Distance = ({moth_ground_speed_m_per_min} / 60 m/s) * {led_propagation_time_s:.1f} s = {additional_distance_m} m\n")

# --- Step 4: Calculate the total displacement ---
# The total displacement is the 1m to the halfway point plus the additional distance.
distance_to_halfway_m = tunnel_length_m / 2.0
total_displacement_m = distance_to_halfway_m + additional_distance_m

print("The moth's total displacement is the 1m to the halfway point plus this additional distance.")
print("Final Equation:")
print(f"Total Displacement = Distance to Halfway + ( (Moth Air Speed - Wind Speed) / 60 * (Number of LEDs - 1) * LED Delay )")
print(f"Total Displacement = {distance_to_halfway_m} m + ( ({moth_air_speed_m_per_min} m/min - {wind_speed_m_per_min} m/min) / 60 s/min * ({num_leds} - 1) * {led_delay_s} s )")
print(f"Total Displacement = {distance_to_halfway_m} m + {additional_distance_m} m = {total_displacement_m} m")
print(f"\nThe moth's displacement is {total_displacement_m} m, which is approximately {total_displacement_m:.2f} m.")
