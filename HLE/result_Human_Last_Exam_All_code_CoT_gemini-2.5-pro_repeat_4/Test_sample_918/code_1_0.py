import math

# Define the constants from the problem statement
tunnel_length_m = 2.0
moth_airspeed_mpm = 5.675  # m/min
wind_speed_mpm = 5.0      # m/min
num_leds = 80
led_delay_s = 0.3         # seconds

# --- Step 1: Calculate the moth's speed relative to the ground ---
# The moth flies west, and the wind blows east, so their effects are opposing.
moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
print("Step 1: Calculate the moth's ground speed.")
print(f"The moth's speed in air ({moth_airspeed_mpm} m/min) is reduced by the headwind ({wind_speed_mpm} m/min).")
print(f"Moth's ground speed = {moth_airspeed_mpm} - {wind_speed_mpm} = {moth_ground_speed_mpm:.3f} m/min\n")

# Convert ground speed to meters per second for consistency with other timings.
moth_ground_speed_mps = moth_ground_speed_mpm / 60.0
print("Converting speed to meters per second:")
print(f"Moth's ground speed = {moth_ground_speed_mpm:.3f} m/min / 60 s/min = {moth_ground_speed_mps:.5f} m/s\n")


# --- Step 2: Calculate the time it takes for the moth to reach the halfway point ---
# The LED blinking sequence starts when the moth is halfway through the tunnel.
halfway_distance_m = tunnel_length_m / 2.0
time_to_halfway_s = halfway_distance_m / moth_ground_speed_mps
print("Step 2: Calculate the time until the LED sequence starts.")
print(f"The sequence starts when the moth reaches the halfway point ({halfway_distance_m} m).")
print(f"Time to reach halfway point = {halfway_distance_m} m / {moth_ground_speed_mps:.5f} m/s = {time_to_halfway_s:.3f} s\n")


# --- Step 3: Calculate the time for the light wave to travel across the tunnel ---
# The blink propagates from the westernmost LED to the easternmost.
# With 80 LEDs, there are 79 intervals between them.
num_intervals = num_leds - 1
time_for_light_wave_s = num_intervals * led_delay_s
print("Step 3: Calculate the time for the light wave to propagate.")
print(f"The wave travels across {num_intervals} LED intervals, with a delay of {led_delay_s} s each.")
print(f"Time for light wave = {num_intervals} * {led_delay_s} s = {time_for_light_wave_s:.1f} s\n")


# --- Step 4: Calculate the total time and final displacement ---
# The total time is the sum of the time to reach the halfway point and the light wave propagation time.
total_time_s = time_to_halfway_s + time_for_light_wave_s
final_displacement_m = moth_ground_speed_mps * total_time_s
print("Step 4: Calculate the final displacement.")
print(f"The moth flies for a total time of {time_to_halfway_s:.3f} s + {time_for_light_wave_s:.1f} s = {total_time_s:.3f} s.")
print("\n--- FINAL CALCULATION ---")
print(f"Displacement = Ground Speed * Total Time")
print(f"Displacement = {moth_ground_speed_mps:.5f} m/s * {total_time_s:.3f} s")
print(f"Final Displacement = {final_displacement_m:.2f} m")