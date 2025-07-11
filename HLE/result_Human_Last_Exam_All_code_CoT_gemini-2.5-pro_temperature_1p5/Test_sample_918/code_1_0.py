import math

# Define the given constants
tunnel_length_m = 2.0
wind_speed_m_per_min = 5.0
moth_air_speed_m_per_min = 5.675
num_lights = 80
light_delay_s = 0.3

# --- Step 1: Calculate the Moth's Effective Speed ---
# The moth flies against the wind, so its speed is reduced.
moth_effective_speed_m_per_min = moth_air_speed_m_per_min - wind_speed_m_per_min

# --- Step 2: Determine the Moth's Position When the Lights Start Blinking ---
# The blinking starts when the moth is halfway through the tunnel.
trigger_displacement_m = tunnel_length_m / 2.0

# --- Step 3: Calculate the Time for the LED Cascade ---
# The cascade travels from the 1st light (west) to the 80th light (east).
# There are num_lights - 1 intervals between the lights.
num_intervals = num_lights - 1
cascade_duration_s = num_intervals * light_delay_s
# Convert the cascade duration from seconds to minutes for consistency.
cascade_duration_min = cascade_duration_s / 60.0

# --- Step 4: Calculate the Distance Traveled by the Moth During the Cascade ---
# The moth continues to fly during the time it takes the light to cascade across the tunnel.
distance_traveled_during_cascade_m = moth_effective_speed_m_per_min * cascade_duration_min

# --- Step 5: Calculate the Final Total Displacement ---
# The final displacement is the initial 1m to trigger the lights plus the distance traveled during the cascade.
final_displacement_m = trigger_displacement_m + distance_traveled_during_cascade_m

print("The moth's final displacement is the distance it traveled to trigger the lights plus the distance it traveled while the light cascade occurred.")
print("\nEquation for the final displacement:")
print(f"Displacement = (Tunnel Length / 2) + (Moth's Effective Speed * Cascade Time)")
# The equation with numbers:
# Displacement = 1.0 m + (0.675 m/min * ((80 - 1) * 0.3 s / 60 s/min))
# To show the values clearly, we break it down
print("\nCalculation with values:")
print(f"Displacement = {trigger_displacement_m:.2f} m + ({moth_effective_speed_m_per_min:.3f} m/min * (({num_lights} - 1) * {light_delay_s:.1f} s / 60 s/min))")
print(f"Displacement = {trigger_displacement_m:.2f} m + ({moth_effective_speed_m_per_min:.3f} m/min * {cascade_duration_min:.4f} min)")
print(f"Displacement = {trigger_displacement_m:.2f} m + {distance_traveled_during_cascade_m:.4f} m")
print(f"Final Displacement = {final_displacement_m:.4f} m")

# Round the answer for the multiple choice options
final_answer_rounded = round(final_displacement_m, 2)
print(f"\nThe moth's displacement when the LED nearest the eastern end of the tunnel blinks is approximately {final_answer_rounded} m.")
<<<D>>>