import math

# Problem parameters
tunnel_length = 2.0  # in meters
moth_start_position_desc = "eastern end"

wind_speed_mpm = 5.0  # in meters per minute
moth_air_speed_mpm = 5.675  # in meters per minute

num_leds = 80
led_blink_delay_s = 0.3  # in seconds

# Step 1: Calculate the moth's ground speed.
# The moth flies against the wind, so we subtract the wind speed from its air speed.
moth_ground_speed_mpm = moth_air_speed_mpm - wind_speed_mpm
# Convert ground speed to meters per second for calculations with time in seconds.
moth_ground_speed_mps = moth_ground_speed_mpm / 60.0

# Step 2: Calculate the time for the light signal to propagate across the LEDs.
# The signal travels across 79 intervals between the 80 LEDs.
num_intervals = num_leds - 1
propagation_time_s = num_intervals * led_blink_delay_s

# Step 3: Calculate the total displacement.
# The total displacement is the distance to the halfway point (1m) plus the
# additional distance the moth flies during the light propagation time.
distance_to_halfway_m = tunnel_length / 2.0
additional_distance_m = moth_ground_speed_mps * propagation_time_s
total_displacement_m = distance_to_halfway_m + additional_distance_m

# Final Output: Print the final equation with all the numbers.
print("The moth's displacement is the distance to the halfway point plus the extra distance traveled during the LED signal propagation.")
print(f"Total Displacement = (Tunnel Length / 2) + ((Moth Air Speed - Wind Speed) / 60) * (Number of LEDs - 1) * (LED Delay)")
print(f"Total Displacement = ({tunnel_length} / 2) + (({moth_air_speed_mpm} - {wind_speed_mpm}) / 60) * ({num_leds} - 1) * {led_blink_delay_s} = {total_displacement_m:.2f}m")

<<<D>>>