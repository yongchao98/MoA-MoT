def calculate_moth_displacement():
    """
    Calculates the moth's displacement based on the problem description.
    """
    # --- Constants from the problem statement ---
    tunnel_length_m = 2.0
    wind_speed_mpm = 5.0  # m/min
    moth_airspeed_mpm = 5.675  # m/min
    num_leds = 80
    led_delay_s = 0.3  # seconds
    seconds_per_minute = 60.0

    # --- Step 1: Calculate the moth's net speed ---
    # The moth flies west against the eastward wind.
    moth_net_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

    # --- Step 2: The blinking starts when the moth is at the halfway point ---
    # The displacement at this point is half the tunnel length.
    distance_to_halfway_m = tunnel_length_m / 2.0

    # --- Step 3: Calculate the time for the LED blink sequence to complete ---
    # The time it takes for the blink to travel from the 1st to the 80th LED.
    # This involves (80 - 1) intervals.
    propagation_time_s = (num_leds - 1) * led_delay_s
    # Convert propagation time to minutes to use with the moth's speed.
    propagation_time_min = propagation_time_s / seconds_per_minute

    # --- Step 4: Calculate the additional distance traveled by the moth ---
    # During the LED propagation, the moth continues to fly.
    additional_distance_m = moth_net_speed_mpm * propagation_time_min

    # --- Step 5: Calculate the total displacement ---
    # Total displacement = distance to halfway + additional distance
    final_displacement_m = distance_to_halfway_m + additional_distance_m

    # --- Final Output ---
    print("The moth's total displacement is the distance to the halfway point plus the additional distance traveled while the LEDs blink across the tunnel.")
    print("Equation: Displacement = (Tunnel Length / 2) + (Moth Airspeed - Wind Speed) * ((Number of LEDs - 1) * Delay between LEDs) / Seconds per Minute")
    print(f"Calculation: {distance_to_halfway_m} + ({moth_airspeed_mpm} - {wind_speed_mpm}) * (({num_leds} - 1) * {led_delay_s}) / {seconds_per_minute}")
    print(f"Result: The moth's final displacement is {final_displacement_m:.6f} meters.")
    print(f"This rounds to {final_displacement_m:.2f}m.")

calculate_moth_displacement()