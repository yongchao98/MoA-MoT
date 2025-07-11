import math

# Define the given parameters
tunnel_length_m = 2.0
wind_speed_mpmin = 5.0  # Eastward
moth_airspeed_mpmin = 5.675
num_leds = 80
led_blink_delay_s = 0.3

# --- Step 1: Calculate the moth's speed relative to the ground ---
# The moth flies west, against the eastward wind.
# We'll consider the westward direction as the direction of displacement.
moth_ground_speed_mpmin = moth_airspeed_mpmin - wind_speed_mpmin

# --- Step 2: Calculate the time to reach the halfway point ---
# The blinking starts when the moth has traveled 1 meter.
distance_to_halfway_m = tunnel_length_m / 2.0
time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpmin

# --- Step 3: Calculate the time for the blink signal to travel ---
# The signal travels across the 79 gaps between the 80 LEDs.
signal_travel_time_s = (num_leds - 1) * led_blink_delay_s
# Convert signal travel time to minutes for consistency
signal_travel_time_min = signal_travel_time_s / 60.0

# --- Step 4: Calculate the total displacement ---
# The total displacement is the distance to the halfway point (1m) plus the
# additional distance the moth travels while the blink signal propagates.
additional_distance_m = moth_ground_speed_mpmin * signal_travel_time_min
total_displacement_m = distance_to_halfway_m + additional_distance_m

# --- Step 5: Print the final equation and result ---
# The final equation is:
# Displacement = (Distance to Halfway) + (Moth's Ground Speed * Signal Travel Time)
print("The final equation for the moth's displacement is:")
print("Displacement = (Distance to Halfway) + (Moth's Ground Speed * Signal Travel Time)")
print(f"Displacement = {distance_to_halfway_m} m + (({moth_airspeed_mpmin} m/min - {wind_speed_mpmin} m/min) * (({num_leds} - 1) * {led_blink_delay_s} s / 60 s/min))")
print(f"Displacement = {distance_to_halfway_m} m + ({moth_ground_speed_mpmin:.3f} m/min * {signal_travel_time_min:.3f} min)")
print(f"Displacement = {distance_to_halfway_m} m + {additional_distance_m:.6f} m")
print(f"Displacement = {total_displacement_m:.6f} m")
print(f"\nThe moth's displacement when the easternmost LED blinks is approximately {total_displacement_m:.2f}m.")
<<<D>>>