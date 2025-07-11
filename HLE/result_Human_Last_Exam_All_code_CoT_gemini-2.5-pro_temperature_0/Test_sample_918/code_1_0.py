import math

# Step 1: Define constants and calculate the moth's ground speed.
tunnel_length_m = 2.0
air_speed_mpm = 5.0
moth_airspeed_mpm = 5.675
num_leds = 80
led_delay_s = 0.3

# Convert speeds from m/min to m/s
air_speed_mps = air_speed_mpm / 60.0
moth_airspeed_mps = moth_airspeed_mpm / 60.0

# The moth flies west (against the wind), so its ground speed is its speed minus the wind speed.
moth_ground_speed_mps = moth_airspeed_mps - air_speed_mps

# Step 2: Calculate the displacement until the blinking starts.
# The blinking starts when the moth is halfway, which is a displacement of 1m.
displacement_at_blink_start = tunnel_length_m / 2.0

# Step 3: Calculate the time it takes for the blink to travel from the west end to the east end.
# This is the time the moth continues to travel after the blinking sequence starts.
time_for_blink_to_travel_s = (num_leds - 1) * led_delay_s

# Step 4: Calculate the additional displacement during the blink travel time.
additional_displacement_m = moth_ground_speed_mps * time_for_blink_to_travel_s

# Step 5: Calculate the total displacement.
# Total Displacement = Displacement until blinking starts + Additional displacement during blink travel
total_displacement_m = displacement_at_blink_start + additional_displacement_m

# Output the explanation and the final equation with numbers.
print("The moth's total displacement is the distance it has traveled from the start.")
print("This can be calculated as the sum of two parts:")
print("1. The distance to the halfway point (where the blinking starts).")
print("2. The additional distance the moth travels while the blink signal propagates from the west to the east end.")
print("\nThe final equation is:")
print("Total Displacement = (Distance to Halfway) + (Moth Ground Speed) * (Blink Travel Time)")
print("\nPlugging in the numbers:")
print(f"Total Displacement = {displacement_at_blink_start} m + ({moth_ground_speed_mps:.5f} m/s) * (({num_leds} - 1) * {led_delay_s} s)")
print(f"Total Displacement = {displacement_at_blink_start} m + ({moth_ground_speed_mps:.5f} m/s) * ({time_for_blink_to_travel_s} s)")
print(f"Total Displacement = {displacement_at_blink_start} m + {additional_displacement_m:.6f} m")
print(f"Total Displacement = {total_displacement_m:.6f} m")
print(f"\nThe moth's displacement is {total_displacement_m:.2f}m.")
