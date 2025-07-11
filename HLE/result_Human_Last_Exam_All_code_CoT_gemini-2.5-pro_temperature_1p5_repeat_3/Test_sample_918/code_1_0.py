import math

# Define the given parameters
tunnel_length_m = 2.0
wind_speed_mpm = 5.0  # meters per minute
moth_airspeed_mpm = 5.675  # meters per minute
num_leds = 80
led_blink_delay_s = 0.3  # seconds
seconds_per_minute = 60.0

# 1. Calculate the moth's speed relative to the ground (ground speed)
# The moth flies against the wind, so their speeds are subtracted.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

# 2. Calculate the time it takes for the moth to reach the halfway point
# The distance to the halfway point is half the tunnel's length.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm

# 3. Calculate the time it takes for the light signal to travel from west to east
# The signal crosses (number of LEDs - 1) gaps.
time_for_light_travel_s = (num_leds - 1) * led_blink_delay_s
time_for_light_travel_min = time_for_light_travel_s / seconds_per_minute

# 4. Calculate the total time the moth has been flying
# This is the time to reach the halfway point plus the time for the light to travel.
total_time_min = time_to_halfway_min + time_for_light_travel_min

# 5. Calculate the total displacement of the moth from its starting point
# Displacement = Speed * Total Time
total_displacement_m = moth_ground_speed_mpm * total_time_min

# Print the final equation with each component value
# The calculation can be simplified to: 1 meter (from the first leg of the journey)
# plus the distance traveled during the light signal's travel time.
# Displacement = 1 + (moth_ground_speed_mpm * time_for_light_travel_min)
print("The moth's displacement is calculated by its ground speed multiplied by the total time of flight.")
print("Final Equation:")
print(f"{moth_ground_speed_mpm:.3f} m/min * (({distance_to_halfway_m:.1f}m / {moth_ground_speed_mpm:.3f} m/min) + ({time_for_light_travel_s:.1f}s / {seconds_per_minute:.1f} s/min)) = {total_displacement_m:.2f} m")
<<<D>>>