import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the problem description.
    """
    # Step 1: Define variables from the problem
    tunnel_length_m = 2.0
    wind_speed_mpm = 5.0  # m/min, flowing east
    moth_airspeed_mpm = 5.675  # m/min
    led_count = 80
    led_delay_s = 0.3  # seconds

    # Convert units for consistency (using minutes for time)
    led_delay_min = led_delay_s / 60.0

    # Step 2: Calculate the moth's speed relative to the ground
    # The moth flies west (against the wind), so its ground speed is its air speed minus the wind speed.
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

    # Step 3: Calculate the time for the moth to reach the halfway point
    # The blinking starts when the moth has traveled half the tunnel's length.
    distance_to_halfway_m = tunnel_length_m / 2.0
    time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm

    # Step 4: Calculate the time for the blink to travel from the west end to the east end
    # The blink signal propagates across (led_count - 1) intervals.
    time_propagation_min = (led_count - 1) * led_delay_min

    # Step 5: Calculate the total time elapsed from the moth's start
    total_time_min = time_to_halfway_min + time_propagation_min

    # Step 6: Calculate the moth's total displacement from its starting point
    total_displacement_m = moth_ground_speed_mpm * total_time_min

    # --- Output the results ---
    print("Step 1: Calculating the moth's speed relative to the ground.")
    print(f"The moth's effective speed is {moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min = {moth_ground_speed_mpm} m/min westward.\n")

    print("Step 2: Calculating the time until the LED sequence begins.")
    print(f"The moth reaches the halfway point (1.0m) in {distance_to_halfway_m} m / {moth_ground_speed_mpm} m/min = {time_to_halfway_min:.4f} minutes.\n")
    
    print("Step 3: Calculating the time for the LED signal to propagate.")
    propagation_time_s = (led_count - 1) * led_delay_s
    print(f"The signal travels across {led_count - 1} LED gaps in ({led_count - 1}) * {led_delay_s}s = {propagation_time_s:.1f} seconds, which is {time_propagation_min:.4f} minutes.\n")

    print("Step 4: Calculating the total displacement.")
    print(f"The total time elapsed from the moth's start is {time_to_halfway_min:.4f} min + {time_propagation_min:.4f} min = {total_time_min:.4f} minutes.")
    
    print("\nThe final equation for the displacement is: ground_speed * (time_to_halfway + time_propagation)")
    print(f"Displacement = {moth_ground_speed_mpm} m/min * ({time_to_halfway_min:.4f} min + {time_propagation_min:.4f} min)")
    print(f"Displacement = {moth_ground_speed_mpm} m/min * {total_time_min:.4f} min")
    print(f"Final Displacement = {total_displacement_m:.4f} m\n")
    
    print(f"The calculated displacement of {total_displacement_m:.4f}m rounds to 1.27m.")

solve_moth_displacement()