import math

# Define the given parameters from the problem description
tunnel_length_m = 2.0
moth_airspeed_m_per_min = 5.675
wind_speed_m_per_min = 5.0
num_leds = 80
time_delay_s = 0.3

# Step 1: Calculate the moth's effective speed relative to the tunnel.
# The moth flies west, and the wind blows east, so their speeds are in opposition.
effective_speed_m_per_min = moth_airspeed_m_per_min - wind_speed_m_per_min
# Convert the speed from meters per minute to meters per second for consistency.
effective_speed_m_per_s = effective_speed_m_per_min / 60.0

# Step 2: Calculate the total time for the LED sequence to propagate.
# The signal travels from the 1st LED to the 80th, crossing 79 intervals.
num_intervals = num_leds - 1
sequence_duration_s = num_intervals * time_delay_s

# Step 3: Calculate the distance the moth travels during the LED sequence.
# The sequence starts when the moth is at the halfway point. We need to find
# how much further it flies from that point.
distance_after_halfway_m = effective_speed_m_per_s * sequence_duration_s

# Step 4: Calculate the total displacement from the starting point (eastern end).
# The displacement is the distance to the halfway point plus the additional distance traveled.
distance_to_halfway_m = tunnel_length_m / 2.0
total_displacement_m = distance_to_halfway_m + distance_after_halfway_m

# Print the explanation and the final equation with all the numbers.
print("The moth's displacement is calculated as the initial distance to the halfway point plus the additional distance it travels while the LED signal propagates from the first to the last light.")
print("\nFinal Equation:")
print(f"Displacement = (Tunnel Length / 2) + ((Moth Speed - Wind Speed) / 60) * ((Number of LEDs - 1) * Time Delay)")
print(f"Displacement = ({tunnel_length_m}m / 2) + (({moth_airspeed_m_per_min}m/min - {wind_speed_m_per_min}m/min) / 60s/min) * (({num_leds} - 1) * {time_delay_s}s)")

print("\nStep-by-step Calculation:")
print(f"Displacement = {distance_to_halfway_m}m + ({effective_speed_m_per_min}m/min / 60s/min) * ({sequence_duration_s}s)")
print(f"Displacement = {distance_to_halfway_m}m + {effective_speed_m_per_s:.5f}m/s * {sequence_duration_s}s")
print(f"Displacement = {distance_to_halfway_m}m + {distance_after_halfway_m:.5f}m")
print(f"Total Displacement = {total_displacement_m:.5f}m")

# Round the final answer to two decimal places to match the options.
final_answer = round(total_displacement_m, 2)
print(f"\nThe moth's displacement when the LED nearest the eastern end blinks is approximately {final_answer}m.")