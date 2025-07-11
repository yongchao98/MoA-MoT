# Constants from the problem description
tunnel_length_m = 2.0
wind_speed_mpm = 5.0
moth_airspeed_mpm = 5.675
num_leds = 80
led_delay_s = 0.3
seconds_per_minute = 60

# Step 1: Calculate the moth's net speed relative to the ground
# The moth flies west against the eastward wind.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

# Step 2: Calculate the time for the moth to reach the halfway point
halfway_distance_m = tunnel_length_m / 2.0

# Using the simplified equation: Displacement = half_distance + speed * light_travel_time
# This avoids intermediate rounding. First, calculate the light travel time.

# Step 3: Calculate the time for the LED signal to travel from west to east
# The signal propagates across (num_leds - 1) intervals.
light_travel_time_s = (num_leds - 1) * led_delay_s
# Convert light travel time to minutes to match the speed unit
light_travel_time_min = light_travel_time_s / seconds_per_minute

# Step 4 & 5: Calculate the final displacement
# The total displacement is the distance to the halfway point (where the light event is triggered)
# plus the additional distance the moth travels during the light signal's propagation.
additional_distance_m = moth_ground_speed_mpm * light_travel_time_min
final_displacement_m = halfway_distance_m + additional_distance_m

# Print the final equation with all the numbers
print("The final displacement is calculated by adding the halfway distance to the additional distance the moth travels while the lights are blinking across the tunnel.")
print("Final Displacement = Halfway Distance + (Moth's Ground Speed * Light Signal Travel Time)")
print(f"Final Displacement = {halfway_distance_m}m + ({moth_ground_speed_mpm} m/min * {light_travel_time_min} min)")
print(f"Final Displacement = {halfway_distance_m}m + {additional_distance_m}m")
print(f"Final Displacement = {final_displacement_m:.4f}m")
print(f"Rounding to two decimal places, the displacement is {final_displacement_m:.2f}m.")
