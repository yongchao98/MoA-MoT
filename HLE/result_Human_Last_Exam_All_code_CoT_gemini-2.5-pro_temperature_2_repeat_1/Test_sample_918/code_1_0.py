import math

def solve_moth_problem():
    """
    Calculates the moth's displacement based on the problem's parameters.
    """
    # Given parameters
    tunnel_length_m = 2.0
    moth_airspeed_mpm = 5.675  # meters per minute
    wind_speed_mpm = 5.0      # meters per minute
    num_leds = 80
    led_blink_delay_s = 0.3   # seconds

    # Step 1: Calculate the moth's ground speed in m/s
    # The moth is flying west, against the eastward wind.
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
    moth_ground_speed_ms = moth_ground_speed_mpm / 60.0

    print("Step 1: Calculating the moth's speed relative to the ground.")
    print(f"Moth's ground speed = {moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min = {moth_ground_speed_mpm} m/min")
    print(f"Moth's ground speed in m/s = {moth_ground_speed_mpm} / 60 = {moth_ground_speed_ms:.5f} m/s\n")

    # Step 2: Calculate the time for the moth to reach the tunnel's midpoint
    # The blinking starts when the moth is halfway, i.e., has traveled 1m.
    distance_to_midpoint_m = tunnel_length_m / 2.0
    time_to_midpoint_s = distance_to_midpoint_m / moth_ground_speed_ms

    print("Step 2: Calculating the time to reach the midpoint.")
    print(f"Time to midpoint = {distance_to_midpoint_m} m / {moth_ground_speed_ms:.5f} m/s = {time_to_midpoint_s:.2f} s\n")

    # Step 3: Calculate the time for the LED blink to propagate from west to east
    # The signal travels from LED 1 to LED 80, so there are 79 gaps.
    signal_propagation_time_s = (num_leds - 1) * led_blink_delay_s

    print("Step 3: Calculating LED signal propagation time.")
    print(f"Signal time = ({num_leds} - 1) * {led_blink_delay_s} s = {signal_propagation_time_s:.1f} s\n")

    # Step 4: Calculate the total displacement (distance traveled)
    # Total distance = distance to midpoint + distance traveled during signal propagation
    distance_after_midpoint_m = moth_ground_speed_ms * signal_propagation_time_s
    total_displacement_m = distance_to_midpoint_m + distance_after_midpoint_m
    
    print("Step 4: Calculating the final displacement.")
    print("The moth's total displacement is the sum of the distance to the midpoint and the distance traveled while the light signal propagates.")
    print("\n--- FINAL EQUATION ---")
    print(f"Displacement = (Distance to Midpoint) + (Moth Ground Speed * Signal Propagation Time)")
    print(f"Displacement = {distance_to_midpoint_m} m + ({moth_ground_speed_ms:.5f} m/s * {signal_propagation_time_s:.1f} s)")
    print(f"Displacement = {distance_to_midpoint_m} m + {distance_after_midpoint_m:.6f} m")
    print(f"Displacement = {total_displacement_m:.6f} m")

    print(f"\nThe moth's displacement when the LED nearest the eastern end of the tunnel blinks is approximately {total_displacement_m:.2f} m.")

solve_moth_problem()