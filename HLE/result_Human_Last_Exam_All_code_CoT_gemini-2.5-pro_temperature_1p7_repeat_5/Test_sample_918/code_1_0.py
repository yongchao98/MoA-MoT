def solve_moth_problem():
    """
    This script calculates the moth's displacement based on the provided conditions.
    """
    # 1. Define the given variables
    tunnel_length_m = 2.0
    wind_speed_mpm = 5.0      # m/min
    moth_airspeed_mpm = 5.675 # m/min
    num_leds = 80
    led_delay_s = 0.3         # seconds

    # 2. Calculate the moth's ground speed
    # The moth flies west, and the wind blows east, so their speeds subtract.
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
    
    # 3. Calculate the distance to the halfway point (trigger point)
    halfway_distance_m = tunnel_length_m / 2.0

    # 4. Calculate the total time for the LED signal to propagate from west to east
    # There are (num_leds - 1) intervals between the 80 LEDs.
    led_intervals = num_leds - 1
    propagation_time_s = led_intervals * led_delay_s

    # 5. Calculate the additional distance the moth travels during the propagation time
    # First, convert moth's ground speed from m/min to m/s
    moth_ground_speed_mps = moth_ground_speed_mpm / 60.0
    # Then, calculate the distance: distance = speed * time
    additional_distance_m = moth_ground_speed_mps * propagation_time_s

    # 6. Calculate the final displacement
    # This is the initial distance to the halfway point plus the additional distance traveled.
    final_displacement_m = halfway_distance_m + additional_distance_m
    
    # 7. Print the breakdown of the final equation and the result
    print("The moth's final displacement is the sum of the distance to the halfway point and the extra distance traveled while the LED signal propagates.")
    print("\n--- Calculation Breakdown ---")
    print(f"Moth's ground speed: {moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min = {moth_ground_speed_mpm} m/min")
    print(f"Distance to halfway point (event trigger): {halfway_distance_m} m")
    print(f"LED signal propagation time: ({num_leds} - 1) * {led_delay_s}s = {propagation_time_s}s")
    print(f"Additional distance traveled by moth: ({moth_ground_speed_mpm} / 60) m/s * {propagation_time_s}s = {additional_distance_m:.5f} m")

    print("\n--- Final Equation ---")
    print(f"Final Displacement = (Distance to halfway) + (Ground Speed * Propagation Time)")
    # We output each number in the final equation as requested
    print(f"Final Displacement = {halfway_distance_m}m + {additional_distance_m:.5f}m")
    print(f"Final Displacement = {final_displacement_m:.5f}m")
    
    # The result 1.26663 is closest to 1.27m
    print("\nThis rounds to 1.27m.")


solve_moth_problem()
<<<D>>>