import math

# Define the initial parameters from the problem
tunnel_length_m = 2.0
wind_speed_m_per_min = 5.0
moth_speed_air_m_per_min = 5.675
num_leds = 80
led_delay_s = 0.3
seconds_per_minute = 60.0

# --- Step 1: Calculate the moth's effective speed ---
# The moth flies west (against the wind), so its speed relative to the ground is its own speed minus the wind speed.
effective_speed_m_per_min = moth_speed_air_m_per_min - wind_speed_m_per_min
print("--- Calculation Steps ---")
print(f"Step 1: Calculate the moth's effective speed relative to the ground.")
print(f"  Moth's speed in still air (westward): {moth_speed_air_m_per_min} m/min")
print(f"  Wind speed (eastward): {wind_speed_m_per_min} m/min")
print(f"  Effective speed = {moth_speed_air_m_per_min} - {wind_speed_m_per_min} = {effective_speed_m_per_min:.3f} m/min (westward)\n")


# --- Step 2: Calculate the time for the moth to reach the halfway point ---
# The light sequence starts when the moth is halfway, which is 1m from the start.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / effective_speed_m_per_min
print(f"Step 2: Calculate the time for the moth to trigger the light sequence.")
print(f"  The sequence starts when the moth reaches the halfway point ({distance_to_halfway_m}m).")
print(f"  Time to reach halfway = Distance / Speed = {distance_to_halfway_m}m / {effective_speed_m_per_min:.3f} m/min = {time_to_halfway_min:.4f} minutes\n")


# --- Step 3: Calculate the total time for the LED cascade ---
# The cascade must travel across 79 gaps between the 80 LEDs.
led_gaps = num_leds - 1
cascade_time_s = led_gaps * led_delay_s
cascade_time_min = cascade_time_s / seconds_per_minute
print(f"Step 3: Calculate the time for the light signal to reach the last LED.")
print(f"  The signal travels across {led_gaps} gaps between {num_leds} LEDs.")
print(f"  Total cascade time = {led_gaps} * {led_delay_s}s = {cascade_time_s:.1f} seconds, or {cascade_time_min:.4f} minutes\n")


# --- Step 4: Calculate the total elapsed time ---
# This is the time from the moth's start until the final LED blinks.
total_time_min = time_to_halfway_min + cascade_time_min
print(f"Step 4: Calculate the total time from the moth's start.")
print(f"  Total time = Time to halfway + Cascade time = {time_to_halfway_min:.4f} min + {cascade_time_min:.4f} min = {total_time_min:.4f} minutes\n")


# --- Step 5: Calculate the moth's final displacement ---
# Displacement = effective speed * total time
total_displacement_m = effective_speed_m_per_min * total_time_min
print(f"Step 5: Calculate the moth's total displacement.")
print(f"  Displacement = Effective Speed * Total Time = {effective_speed_m_per_min:.3f} m/min * {total_time_min:.4f} min = {total_displacement_m:.4f} meters")
print(f"  Rounding to two decimal places, the displacement is {total_displacement_m:.2f}m.\n")


# --- Final Equation ---
print("--- Final Equation with Initial Values ---")
print("Displacement = (Moth Speed - Wind Speed) * ( (Tunnel Length / 2) / (Moth Speed - Wind Speed) + ( (Num LEDs - 1) * LED Delay) / 60 )")
print(f"Displacement = ({moth_speed_air_m_per_min} - {wind_speed_m_per_min}) * ( ({tunnel_length_m} / 2) / ({moth_speed_air_m_per_min} - {wind_speed_m_per_min}) + ( ({num_leds} - 1) * {led_delay_s}) / {seconds_per_minute} ) = {total_displacement_m:.5f} m")
