def solve_moth_problem():
    """
    Calculates the moth's displacement based on the problem description.
    """
    # Define variables from the problem
    tunnel_length_m = 2.0
    wind_speed_mpm = 5.0  # meters per minute
    moth_airspeed_mpm = 5.675  # meters per minute
    num_leds = 80
    led_delay_s = 0.3  # seconds

    # --- Step 1: Calculate the moth's ground speed ---
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
    # Convert speed to meters per second for calculations involving seconds
    moth_ground_speed_mps = moth_ground_speed_mpm / 60.0

    # --- Step 2: Calculate the time for the moth to reach the halfway point ---
    distance_to_halfway_m = tunnel_length_m / 2.0
    time_to_halfway_s = distance_to_halfway_m / moth_ground_speed_mps

    # --- Step 3: Calculate the LED signal propagation time ---
    # The signal travels across 79 gaps (between 80 LEDs)
    led_propagation_time_s = (num_leds - 1) * led_delay_s

    # --- Step 4: Calculate the total time elapsed ---
    total_time_s = time_to_halfway_s + led_propagation_time_s

    # --- Step 5: Calculate the final displacement of the moth ---
    final_displacement_m = moth_ground_speed_mps * total_time_s

    # --- Print the results and the final equation with all numbers ---
    print("Intermediate Calculations:")
    print(f"Moth's ground speed = ({moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min) / 60 = {moth_ground_speed_mps:.5f} m/s")
    print(f"Time to reach halfway = {distance_to_halfway_m} m / {moth_ground_speed_mps:.5f} m/s = {time_to_halfway_s:.2f} s")
    print(f"LED propagation time = ({num_leds} - 1) * {led_delay_s} s = {led_propagation_time_s:.1f} s")
    print(f"Total time elapsed = {time_to_halfway_s:.2f} s + {led_propagation_time_s:.1f} s = {total_time_s:.2f} s")
    print("\nFinal Displacement Equation:")
    print(f"Displacement = Speed * Time")
    print(f"Displacement = {moth_ground_speed_mps:.5f} m/s * {total_time_s:.2f} s")
    print(f"Final calculated displacement: {final_displacement_m:.3f} m")
    print("\nThis rounds to 1.27m, which corresponds to answer choice D.")

solve_moth_problem()
<<<D>>>