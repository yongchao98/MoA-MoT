# 1. Define the given constants
moth_airspeed_west = 5.675  # m/min
wind_speed_east = 5.0  # m/min
tunnel_length = 2.0  # m
num_lights = 80
blink_delay_per_light = 0.3  # seconds

# 2. Calculate the moth's effective ground speed
moth_ground_speed_mpm = moth_airspeed_west - wind_speed_east

# The LED blinking sequence starts when the moth is at the halfway point.
displacement_at_blink_start = tunnel_length / 2

# 3. Calculate the time it takes for the blink to travel from the west end to the east end.
# There are (num_lights - 1) gaps between the first and last light.
num_gaps = num_lights - 1
blink_propagation_time_sec = num_gaps * blink_delay_per_light
# Convert this time to minutes to be consistent with the speed unit.
blink_propagation_time_min = blink_propagation_time_sec / 60.0

# 4. Calculate the distance the moth travels during this propagation time.
distance_traveled_during_propagation = moth_ground_speed_mpm * blink_propagation_time_min

# 5. The final displacement is the position when the blinking started plus the extra distance traveled.
final_displacement = displacement_at_blink_start + distance_traveled_during_propagation

# Print out the step-by-step calculation
print("Step 1: Calculate the moth's ground speed.")
print(f"Moth's ground speed = Moth's airspeed - Wind speed = {moth_airspeed_west} m/min - {wind_speed_east} m/min = {moth_ground_speed_mpm} m/min")
print("-" * 30)

print("Step 2: Calculate the time for the light signal to travel the tunnel length.")
print(f"Signal propagation time = (Number of lights - 1) * Delay per light = ({num_lights} - 1) * {blink_delay_per_light} s = {blink_propagation_time_sec} s")
print("-" * 30)

print("Step 3: Calculate the distance the moth travels during the signal propagation.")
print(f"Distance = Moth's ground speed * Signal propagation time = {moth_ground_speed_mpm} m/min * ({blink_propagation_time_sec} / 60) min = {distance_traveled_during_propagation:.5f} m")
print("-" * 30)

print("Step 4: Calculate the moth's final displacement from the starting (eastern) end.")
print("The blinking starts when the moth is at the halfway point (1m).")
print("Final Displacement = Displacement at blink start + Distance traveled during propagation")
print(f"Final Displacement = {displacement_at_blink_start} m + {distance_traveled_during_propagation:.5f} m = {final_displacement:.5f} m")
print("-" * 30)

print(f"The moth's final displacement is {final_displacement:.2f} m.")
print("\nFinal Equation:")
print(f"Displacement = {displacement_at_blink_start} + ({moth_airspeed_west} - {wind_speed_east}) * (({num_lights} - 1) * {blink_delay_per_light} / 60)")
print(f"{final_displacement:.2f} = {displacement_at_blink_start} + ({moth_ground_speed_mpm}) * (({num_gaps}) * {blink_delay_per_light} / 60)")