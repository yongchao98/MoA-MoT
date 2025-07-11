import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the problem's conditions.
    """
    # --- Define constants from the problem statement ---
    tunnel_length_m = 2.0
    air_speed_mpm = 5.0  # m/min (eastward)
    moth_airspeed_mpm = 5.675  # m/min
    num_leds = 80
    led_delay_s = 0.3  # seconds

    # --- Step 1: Calculate the moth's ground speed ---
    # The moth flies west, against the eastward wind.
    moth_ground_speed_mpm = moth_airspeed_mpm - air_speed_mpm

    # --- Step 2: Calculate the time for the moth to reach the halfway point ---
    # The blinking starts when the moth is at the 1m mark.
    halfway_distance_m = tunnel_length_m / 2.0

    # --- Step 3: Calculate the time for the LED wave to travel ---
    # The wave travels across 79 intervals (80 LEDs - 1).
    time_for_wave_s = (num_leds - 1) * led_delay_s
    # Convert wave time to minutes for consistent units.
    time_for_wave_min = time_for_wave_s / 60.0

    # --- Step 4 & 5: Calculate the final displacement ---
    # The total displacement is the distance to the halfway point plus the
    # additional distance the moth travels while the light wave propagates.
    # Displacement = (Distance to halfway) + (Distance traveled during wave)
    # Displacement = halfway_distance_m + (moth_ground_speed_mpm * time_for_wave_min)
    final_displacement_m = halfway_distance_m + (moth_ground_speed_mpm * time_for_wave_min)

    # --- Print the final equation and result ---
    print("The final displacement is the sum of the distance to the halfway point and the distance the moth travels while the light wave propagates.")
    print("\nFinal Equation:")
    print("Displacement = (Tunnel Length / 2) + (Moth Ground Speed * Time for Light Wave)")
    print("\nWith the given numbers:")
    print(f"Displacement = ({tunnel_length_m} / 2) + (({moth_airspeed_mpm} - {air_speed_mpm}) * (({num_leds} - 1) * {led_delay_s} / 60))")
    print(f"Displacement = {halfway_distance_m}m + ({moth_ground_speed_mpm} m/min * {time_for_wave_min:.4f} min)")
    print(f"Displacement = {halfway_distance_m}m + {moth_ground_speed_mpm * time_for_wave_min:.4f}m")
    print(f"\nFinal Displacement = {final_displacement_m:.2f}m")

solve_moth_displacement()
<<<D>>>