import math

# --- Problem Parameters ---
tunnel_length_m = 2.0
moth_speed_mpm = 5.675  # Moth's airspeed in meters per minute
wind_speed_mpm = 5.0    # Wind speed in meters per minute (eastward)
num_leds = 80
led_delay_s = 0.3       # Delay in seconds between consecutive LEDs

# --- Calculations ---

# Step 1: Calculate the moth's effective speed.
# The moth flies west, while the wind blows east, so the speeds are opposing.
effective_speed_mpm = moth_speed_mpm - wind_speed_mpm
# Convert speed to meters per second for consistency.
effective_speed_mps = effective_speed_mpm / 60.0

# Step 2: Calculate the time for the moth to reach the halfway point.
# The moth starts at the eastern end and flies 1m to reach the middle.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_s = distance_to_halfway_m / effective_speed_mps

# Step 3: Calculate the time for the light signal to propagate.
# There are 80 lights, which means there are 79 intervals between them.
num_intervals = num_leds - 1
signal_travel_time_s = num_intervals * led_delay_s

# Step 4: Calculate the total time elapsed.
# This is the sum of the time to get to the halfway point and the signal propagation time.
total_time_s = time_to_halfway_s + signal_travel_time_s

# Step 5: Calculate the moth's total displacement from its starting point.
moth_displacement_m = effective_speed_mps * total_time_s

# --- Output the results ---

print("Step 1: Moth's Effective Speed")
print(f"The moth's effective speed against the wind is {moth_speed_mpm} m/min - {wind_speed_mpm} m/min = {effective_speed_mpm:.3f} m/min.")
print(f"In meters per second, this is {effective_speed_mpm:.3f} / 60 = {effective_speed_mps:.5f} m/s.")
print("-" * 40)

print("Step 2: Time to Reach Halfway Point")
print(f"The LED sequence starts when the moth has traveled {distance_to_halfway_m} m.")
print(f"Time to reach halfway = {distance_to_halfway_m} m / {effective_speed_mps:.5f} m/s = {time_to_halfway_s:.2f} s.")
print("-" * 40)

print("Step 3: LED Signal Propagation Time")
print(f"The signal travels across {num_intervals} intervals, each taking {led_delay_s} s.")
print(f"Total signal time = {num_intervals} * {led_delay_s} s = {signal_travel_time_s:.1f} s.")
print("-" * 40)

print("Step 4: Total Time Elapsed")
print(f"Total time until eastern LED blinks = {time_to_halfway_s:.2f} s + {signal_travel_time_s:.1f} s = {total_time_s:.2f} s.")
print("-" * 40)

print("Step 5: Final Moth Displacement")
print("The moth's displacement is its total distance traveled from the start.")
print("Final Equation: Displacement = Effective Speed * Total Time")
print(f"Displacement = {effective_speed_mps:.5f} m/s * {total_time_s:.2f} s")
print(f"\nThe moth's displacement is: {moth_displacement_m:.2f} m")
