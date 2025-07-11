def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the given parameters.
    """
    # Step 0: Define the constants from the problem description
    tunnel_length_m = 2.0
    air_speed_mpm = 5.0  # m/min
    moth_airspeed_mpm = 5.675  # m/min
    num_leds = 80
    led_delay_s = 0.3  # seconds
    secs_per_min = 60

    print("--- Step 1: Calculate the moth's ground speed ---")
    # The moth flies against the wind, so its ground speed is the difference.
    moth_ground_speed_mpm = moth_airspeed_mpm - air_speed_mpm
    print(f"The moth's speed relative to the ground is its airspeed minus the headwind speed.")
    print(f"Ground Speed = {moth_airspeed_mpm} m/min - {air_speed_mpm} m/min = {moth_ground_speed_mpm:.3f} m/min")
    print("-" * 50)

    print("--- Step 2: Calculate the time to reach the halfway point ---")
    # The blinking starts when the moth is halfway, which is 1m from the start.
    distance_to_halfway_m = tunnel_length_m / 2
    time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm
    time_to_halfway_s = time_to_halfway_min * secs_per_min
    print("The blinking starts when the moth has traveled 1m.")
    print(f"Time to halfway point = {distance_to_halfway_m} m / {moth_ground_speed_mpm:.3f} m/min = {time_to_halfway_min:.4f} min")
    print(f"This is equal to {time_to_halfway_min:.4f} min * {secs_per_min} s/min = {time_to_halfway_s:.4f} s")
    print("-" * 50)

    print("--- Step 3: Calculate the light signal propagation time ---")
    # The signal travels across 79 intervals between 80 LEDs.
    num_intervals = num_leds - 1
    propagation_time_s = num_intervals * led_delay_s
    print("The blink signal travels from the 1st to the 80th LED.")
    print(f"Total Propagation Time = ({num_leds} - 1) intervals * {led_delay_s} s/interval")
    print(f"Total Propagation Time = {num_intervals} * {led_delay_s} = {propagation_time_s:.1f} s")
    print("-" * 50)

    print("--- Step 4: Calculate total elapsed time ---")
    # Total time is the sum of time to halfway and the propagation time.
    total_time_s = time_to_halfway_s + propagation_time_s
    print("The total time is the time for the moth to reach the halfway point plus the signal propagation time.")
    print(f"Total Time = {time_to_halfway_s:.4f} s + {propagation_time_s:.1f} s = {total_time_s:.4f} s")
    print("-" * 50)

    print("--- Step 5: Calculate the moth's final displacement ---")
    # Displacement = ground speed * total time. Convert speed to m/s.
    moth_ground_speed_mps = moth_ground_speed_mpm / secs_per_min
    final_displacement_m = moth_ground_speed_mps * total_time_s
    print("The moth's final displacement is its ground speed multiplied by the total elapsed time.")
    print(f"Moth Ground Speed in m/s = {moth_ground_speed_mpm:.3f} m/min / {secs_per_min} s/min = {moth_ground_speed_mps:.5f} m/s")
    print(f"Final Displacement = {moth_ground_speed_mps:.5f} m/s * {total_time_s:.4f} s = {final_displacement_m:.4f} m")
    print("-" * 50)
    
    print(f"\nThe final displacement of the moth is {final_displacement_m:.2f} m.")

solve_moth_displacement()