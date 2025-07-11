def solve_moth_problem():
    """
    Calculates the moth's displacement based on the provided scenario.
    """

    # --- 1. Define Constants (using m/s for consistency) ---
    tunnel_length_m = 2.0
    
    # Speeds
    air_speed_mps = 5.0 / 60.0  # m/s, moving east (+)
    moth_airspeed_mps = 5.675 / 60.0  # m/s
    
    # LED properties
    num_leds = 80
    led_blink_delay_s = 0.3
    
    # Positions
    moth_start_pos_m = 2.0  # Eastern end
    halfway_pos_m = 1.0

    # --- 2. Calculate Moth's Ground Speed ---
    # The moth flies west (-), while the air moves east (+)
    moth_ground_speed_mps = air_speed_mps - moth_airspeed_mps
    
    print(f"Moth's speed relative to the ground: {abs(moth_ground_speed_mps):.5f} m/s towards the west.")

    # --- 3. Calculate Time to Reach Halfway Point ---
    distance_to_halfway_m = moth_start_pos_m - halfway_pos_m
    time_to_halfway_s = distance_to_halfway_m / abs(moth_ground_speed_mps)

    print(f"Time for moth to travel {distance_to_halfway_m}m to the halfway point: {time_to_halfway_s:.2f} s.")

    # --- 4. Calculate Blink Travel Time ---
    # The blink travels from LED #1 to LED #80, covering 79 intervals
    blink_travel_time_s = (num_leds - 1) * led_blink_delay_s

    print(f"Time for LED signal to travel the tunnel length: {blink_travel_time_s:.2f} s.")

    # --- 5. Calculate Total Time ---
    # Total time is from the moth's start until the east LED blinks
    total_time_s = time_to_halfway_s + blink_travel_time_s

    print(f"Total time elapsed from moth's start: {total_time_s:.2f} s.")
    
    # --- 6. Calculate Final Displacement ---
    # Displacement = speed * time. The magnitude is the distance traveled.
    displacement_m = abs(moth_ground_speed_mps) * total_time_s
    
    print("\n--- Final Calculation ---")
    print(f"Moth's Displacement = Moth Ground Speed * Total Time")
    print(f"Displacement = {abs(moth_ground_speed_mps):.5f} m/s * {total_time_s:.2f} s")
    print(f"Displacement = {displacement_m:.2f} m")

solve_moth_problem()
<<<D>>>