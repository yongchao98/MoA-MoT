# Constants from the problem description
tunnel_length_m = 2.0
moth_speed_air_mpm = 5.675  # m/min
wind_speed_mpm = 5.0        # m/min, in the opposite direction
num_lights = 80
delay_between_lights_s = 0.3

# --- Step 1: Calculate the moth's net speed relative to the ground ---
# The moth flies against the wind, so we subtract the wind speed.
moth_net_speed_mpm = moth_speed_air_mpm - wind_speed_mpm
# Convert speed to meters per second for consistency
moth_net_speed_mps = moth_net_speed_mpm / 60.0

# --- Step 2: Calculate the time for the blink to propagate across the tunnel ---
# The blink travels from the 1st to the 80th light, so there are 79 intervals.
num_intervals = num_lights - 1
blink_propagation_time_s = num_intervals * delay_between_lights_s

# --- Step 3: Calculate the final displacement ---
# The blinking starts when the moth is at the halfway point.
halfway_distance_m = tunnel_length_m / 2.0
# The additional distance the moth travels is its net speed multiplied by the propagation time.
additional_distance_m = moth_net_speed_mps * blink_propagation_time_s
# The final displacement is the halfway distance plus the additional distance.
final_displacement_m = halfway_distance_m + additional_distance_m

# --- Step 4: Print the final equation and result ---
# The final equation is: Displacement = Halfway_Distance + (Net_Speed * Propagation_Time)
print("The moth's net speed is {:.4f} m/s.".format(moth_net_speed_mps))
print("The time for the blink to travel from the west to east end is {:.1f} s.".format(blink_propagation_time_s))
print("\nThe final displacement is calculated as:")
print("Displacement = Halfway Distance + (Net Speed * Blink Propagation Time)")
print("Displacement = {}m + ({:.4f} m/s * {:.1f} s)".format(halfway_distance_m, moth_net_speed_mps, blink_propagation_time_s))
print("Displacement = {}m + {:.5f} m".format(halfway_distance_m, additional_distance_m))
print("Final Displacement = {:.5f} m".format(final_displacement_m))