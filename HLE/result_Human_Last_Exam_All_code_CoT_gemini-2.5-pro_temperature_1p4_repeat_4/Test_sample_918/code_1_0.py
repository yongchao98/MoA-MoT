import math

# Define the given parameters
tunnel_length = 2.0  # meters
wind_speed = 5.0  # m/min
moth_airspeed = 5.675  # m/min
num_leds = 80
led_delay_s = 0.3  # seconds

# --- Step 1: Calculate the moth's net speed relative to the ground ---
# The moth flies west, against the eastbound wind.
moth_ground_speed = moth_airspeed - wind_speed
print(f"Step 1: Calculate the moth's net speed.")
print(f"The moth's speed relative to the ground is its airspeed minus the wind speed.")
print(f"Moth's ground speed = {moth_airspeed} m/min - {wind_speed} m/min = {moth_ground_speed:.3f} m/min")
print("-" * 30)

# --- Step 2: Calculate the time for the moth to reach the halfway point ---
# The moth starts at the 2m mark (eastern end) and the halfway point is at the 1m mark.
distance_to_halfway = tunnel_length / 2.0  # meters
time_to_halfway_min = distance_to_halfway / moth_ground_speed
print(f"Step 2: Calculate the time for the moth to reach the halfway point.")
print(f"The blinking starts when the moth travels {distance_to_halfway} m.")
print(f"Time to halfway point = {distance_to_halfway} m / {moth_ground_speed:.3f} m/min = {time_to_halfway_min:.4f} min")
print("-" * 30)

# --- Step 3: Calculate the LED signal propagation time ---
# The signal travels across 80 LEDs, which means there are 79 intervals.
num_intervals = num_leds - 1
propagation_time_s = num_intervals * led_delay_s
# Convert propagation time to minutes to match the moth's speed unit.
propagation_time_min = propagation_time_s / 60.0
print(f"Step 3: Calculate the time for the LED signal to travel the tunnel length.")
print(f"The signal travels across {num_intervals} intervals, each taking {led_delay_s} s.")
print(f"Propagation time = ({num_leds} - 1) * {led_delay_s} s = {propagation_time_s:.1f} s")
print(f"In minutes, this is {propagation_time_s:.1f} s / 60 s/min = {propagation_time_min:.4f} min")
print("-" * 30)

# --- Step 4: Calculate the total time elapsed ---
# This is the sum of the time to reach the halfway point and the signal propagation time.
total_time_min = time_to_halfway_min + propagation_time_min
print(f"Step 4: Calculate the total time from the moth's start.")
print(f"Total time = Time to halfway + Propagation time")
print(f"Total time = {time_to_halfway_min:.4f} min + {propagation_time_min:.4f} min = {total_time_min:.4f} min")
print("-" * 30)

# --- Step 5: Calculate the moth's final displacement ---
# Displacement is the total distance traveled from the start.
displacement = moth_ground_speed * total_time_min
print(f"Step 5: Calculate the moth's total displacement.")
print(f"Displacement = Moth's ground speed * Total time")
print(f"Displacement = {moth_ground_speed:.3f} m/min * {total_time_min:.4f} min = {displacement:.4f} m")
print("-" * 30)

print(f"The final displacement of the moth is {displacement:.2f}m.")
print("\nFinal Equation:")
print(f"Displacement = ({moth_airspeed} - {wind_speed}) * (({distance_to_halfway} / ({moth_airspeed} - {wind_speed})) + (({num_leds} - 1) * {led_delay_s} / 60))")
print(f"{displacement:.2f} = ({moth_airspeed - wind_speed:.3f}) * (({distance_to_halfway} / {moth_airspeed - wind_speed:.3f}) + ({num_intervals * led_delay_s} / 60))")

