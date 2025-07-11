def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the given parameters.
    """
    # Step 1: Define the given variables
    moth_speed_air_mpm = 5.675  # m/min
    wind_speed_mpm = 5.0        # m/min
    tunnel_length_m = 2.0
    num_leds = 80
    led_delay_s = 0.3           # seconds

    # Step 2: Calculate the moth's ground speed
    # The moth flies west, the wind blows east, so they are in opposite directions.
    moth_ground_speed_mpm = moth_speed_air_mpm - wind_speed_mpm
    # Convert ground speed to m/s for calculations with time in seconds
    moth_ground_speed_mps = moth_ground_speed_mpm / 60

    # Step 3: Calculate the time for the LED signal to travel from west to east
    # The time delay between the 1st and 80th LED is (80 - 1) intervals.
    time_led_travel_s = (num_leds - 1) * led_delay_s

    # Step 4: Calculate the distance the moth travels *after* reaching the halfway point.
    # The LED event starts when the moth is at the halfway point (1m displacement).
    # We need to find how much further the moth travels during the LED propagation time.
    distance_after_halfway_m = moth_ground_speed_mps * time_led_travel_s

    # Step 5: Calculate the total displacement from the start (eastern end).
    # This is the 1m to the halfway point plus the additional distance traveled.
    halfway_point_m = tunnel_length_m / 2
    total_displacement_m = halfway_point_m + distance_after_halfway_m

    # Print the step-by-step calculation
    print("Step 1: Calculate the moth's ground speed (speed relative to the tunnel).")
    print(f"Moth's speed in still air: {moth_speed_air_mpm} m/min")
    print(f"Wind speed (opposite direction): {wind_speed_mpm} m/min")
    print(f"Moth's ground speed = {moth_speed_air_mpm} - {wind_speed_mpm} = {moth_ground_speed_mpm:.3f} m/min")
    print(f"Converting to m/s: {moth_ground_speed_mpm:.3f} / 60 = {moth_ground_speed_mps:.5f} m/s\n")

    print("Step 2: Calculate the time for the LED signal to propagate from west to east.")
    print(f"Time = (Number of LEDs - 1) * Delay per LED")
    print(f"Time = ({num_leds} - 1) * {led_delay_s} s = {time_led_travel_s:.1f} s\n")

    print("Step 3: Calculate the total displacement.")
    print("The LED sequence starts when the moth is at the halfway point (1m).")
    print("We calculate the additional distance the moth travels during the LED propagation time.")
    print(f"Additional Distance = Ground Speed * Time = {moth_ground_speed_mps:.5f} m/s * {time_led_travel_s:.1f} s = {distance_after_halfway_m:.6f} m\n")

    print("Final Calculation:")
    print(f"Total Displacement = (Distance to Halfway Point) + (Additional Distance)")
    print(f"Total Displacement = {halfway_point_m} m + {distance_after_halfway_m:.6f} m = {total_displacement_m:.6f} m")
    print(f"\nThe calculated displacement is approximately {total_displacement_m:.2f}m, which corresponds to answer choice D.")


solve_moth_displacement()
<<<D>>>