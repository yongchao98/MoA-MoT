import math

# Define the constants from the problem
tunnel_length_m = 2.0
moth_speed_air_mpm = 5.675  # m/min
wind_speed_mpm = 5.0  # m/min
led_count = 80
led_delay_s = 0.3  # seconds

# --- Step 1: Calculate the moth's effective speed ---
# Convert speeds from meters per minute to meters per second
m_per_min_to_m_per_s = 1.0 / 60.0
moth_speed_air_mps = moth_speed_air_mpm * m_per_min_to_m_per_s
wind_speed_mps = wind_speed_mpm * m_per_min_to_m_per_s

# Moth flies against the wind, so subtract wind speed
moth_effective_speed_mps = moth_speed_air_mps - wind_speed_mps
print(f"1. Moth's effective speed relative to the ground:")
print(f"   ({moth_speed_air_mpm} m/min - {wind_speed_mpm} m/min) / 60 s/min = {moth_effective_speed_mps:.5f} m/s")
print("-" * 30)

# --- Step 2: Calculate time to reach the halfway point ---
# The blinking trigger is when the moth is halfway, which is 1m from the start
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_s = distance_to_halfway_m / moth_effective_speed_mps
print(f"2. Time for the moth to reach the halfway point (1m):")
print(f"   Time = Distance / Speed = {distance_to_halfway_m} m / {moth_effective_speed_mps:.5f} m/s = {time_to_halfway_s:.2f} s")
print("-" * 30)

# --- Step 3: Calculate total time until the final LED blinks ---
# The blinking propagates from the 1st to the 80th LED, so there are 79 delays.
total_blink_delay_s = (led_count - 1) * led_delay_s
# The total time is the time to reach halfway plus the propagation delay
final_blink_time_s = time_to_halfway_s + total_blink_delay_s
print(f"3. Time until the final (80th) LED blinks:")
print(f"   Blink propagation delay = ({led_count} - 1) * {led_delay_s} s = {total_blink_delay_s:.2f} s")
print(f"   Total time from start = {time_to_halfway_s:.2f} s + {total_blink_delay_s:.2f} s = {final_blink_time_s:.2f} s")
print("-" * 30)

# --- Step 4: Calculate the moth's final displacement ---
# Displacement = effective speed * total time
total_displacement_m = moth_effective_speed_mps * final_blink_time_s
print(f"4. Moth's total displacement at the moment of the final blink:")
print(f"   Displacement = Speed * Time")
print(f"   Displacement = {moth_effective_speed_mps:.5f} m/s * {final_blink_time_s:.2f} s")
print(f"   Displacement = {total_displacement_m:.2f} m")
