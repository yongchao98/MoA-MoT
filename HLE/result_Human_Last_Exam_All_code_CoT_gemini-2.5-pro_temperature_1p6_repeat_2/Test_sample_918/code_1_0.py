import math

# --- Givens ---
tunnel_length_m = 2.0
moth_airspeed_mpm = 5.675  # meters per minute
wind_speed_mpm = 5.0       # meters per minute
num_leds = 80
led_delay_s = 0.3          # seconds

# --- Plan ---
# 1. Calculate the moth's speed relative to the ground.
# The moth flies against the wind.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

# 2. Determine the moth's position when the blinking starts.
# Blinking starts when the moth is halfway.
position_at_blink_start_m = tunnel_length_m / 2.0

# 3. Calculate the time it takes for the blink signal to travel from the west end to the east end.
# There are (num_leds - 1) gaps for the signal to cross.
# First, convert the delay from seconds to minutes.
led_delay_min = led_delay_s / 60.0
# Then, calculate the total travel time.
blink_travel_time_min = (num_leds - 1) * led_delay_min

# 4. Calculate the additional distance the moth flies during this travel time.
additional_distance_m = moth_ground_speed_mpm * blink_travel_time_min

# 5. Calculate the final displacement of the moth.
# It's the sum of its position when the blinks started plus the additional distance it flew.
final_displacement_m = position_at_blink_start_m + additional_distance_m

# --- Output ---
# The final equation is: Displacement = Start Position + (Ground Speed * Blink Travel Time)
# Let's show the full calculation with the numbers.
print("The moth's final displacement is calculated as its position when the blinking starts, plus the additional distance it flies while the blink signal travels across the tunnel.")
print("\nEquation: Final Displacement = Position at Blink Start + (Moth Ground Speed * Blink Travel Time)\n")

print("Substituting the values:")
# We format the string to show the full calculation.
# Using the raw numbers that build up to the final result.
print(f"{position_at_blink_start_m} m + ({moth_ground_speed_mpm} m/min * (({num_leds} - 1) * {led_delay_s} s))")
# To make it clearer with consistent units (minutes):
print(f"= {position_at_blink_start_m} m + ({moth_ground_speed_mpm} m/min * {blink_travel_time_min:.3f} min)")
# Final calculation and result
print(f"= {position_at_blink_start_m} m + {additional_distance_m:.4f} m")
print(f"= {final_displacement_m:.4f} m")

# Rounded to two decimal places for the final answer.
print(f"\nThe moth's displacement when the LED nearest the eastern end blinks is {final_displacement_m:.2f}m.")
