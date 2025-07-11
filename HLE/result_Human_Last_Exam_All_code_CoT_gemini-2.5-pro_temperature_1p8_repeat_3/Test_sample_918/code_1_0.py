import math

def solve_moth_problem():
    """
    Solves the moth displacement problem based on the provided parameters.
    """

    # --- Step 1: Define variables and calculate the moth's ground speed ---
    # Speeds are converted from m/min to m/s by dividing by 60
    moth_airspeed_mps = 5.675 / 60  # Moth's speed in still air (m/s)
    wind_speed_mps = 5.0 / 60      # Wind speed (m/s)

    # The moth flies against the wind, so their speeds subtract.
    # We define the moth's direction of travel (west) as positive.
    moth_ground_speed_mps = moth_airspeed_mps - wind_speed_mps

    # --- Step 2: Define tunnel and LED parameters ---
    tunnel_length_m = 2.0
    num_lights = 80
    led_delay_s = 0.3
    
    # --- Step 3: Calculate the displacement ---
    # The problem can be solved by calculating the total distance flown.
    # The moth's displacement is its total distance traveled from the start.
    # Total distance = (distance to halfway) + (distance traveled during blink propagation)

    # The blinking starts when the moth has traveled 1m (half the tunnel length)
    distance_to_halfway_m = tunnel_length_m / 2

    # The blink signal has to travel across 79 gaps (from light 1 to 80)
    num_gaps = num_lights - 1

    # Distance traveled during the blink propagation is (speed * time)
    distance_after_halfway_m = moth_ground_speed_mps * num_gaps * led_delay_s
    
    # Total displacement is the sum of the two distances
    total_displacement_m = distance_to_halfway_m + distance_after_halfway_m

    # --- Step 4: Print the equation and the final result ---
    print("The moth's displacement is the total distance it has flown from the start.")
    print("Displacement = (Distance to halfway) + (Distance flown while LEDs blink to the end)")
    print(f"Displacement = {distance_to_halfway_m:.2f} m + (Moth Ground Speed * (Number of Lights - 1) * LED Delay)")
    print(f"Displacement = {distance_to_halfway_m:.2f} m + ({moth_ground_speed_mps:.5f} m/s * ({num_lights} - 1) * {led_delay_s} s)")
    print(f"Displacement = {distance_to_halfway_m:.2f} m + {distance_after_halfway_m:.5f} m")
    print(f"Final Displacement = {total_displacement_m:.5f} m")
    print(f"Rounded to two decimal places, the displacement is {total_displacement_m:.2f}m.")


solve_moth_problem()
<<<D>>>