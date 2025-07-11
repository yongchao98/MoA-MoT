import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the given parameters.
    """
    # Define the problem parameters from the user's request
    tunnel_length_m = 2.0
    air_speed_m_per_min = 5.0  # Moving east
    moth_airspeed_m_per_min = 5.675 # Moth's speed relative to air, flying west
    num_leds = 80
    led_delay_s = 0.3
    seconds_per_minute = 60.0

    # Step 1: Calculate the moth's speed relative to the ground.
    # The moth flies west, against the eastward wind.
    moth_ground_speed_m_per_min = moth_airspeed_m_per_min - air_speed_m_per_min

    # Step 2: The displacement to trigger the lights is half the tunnel length.
    distance_to_trigger_m = tunnel_length_m / 2.0

    # Step 3: Calculate the time it takes for the light pulse to travel the length of the tunnel.
    # There are (num_leds - 1) intervals between the first and last LED.
    pulse_travel_time_s = (num_leds - 1) * led_delay_s
    pulse_travel_time_min = pulse_travel_time_s / seconds_per_minute

    # Step 4: Calculate the additional distance the moth travels during the light pulse propagation.
    additional_displacement_m = moth_ground_speed_m_per_min * pulse_travel_time_min

    # The total displacement is the distance to trigger the lights plus the additional displacement.
    total_displacement_m = distance_to_trigger_m + additional_displacement_m

    # Print the final equation and the result, showing all the numbers involved.
    print(f"The moth's displacement is the 1m traveled to trigger the lights, plus the distance it travels while the light pulse propagates.")
    print(f"Final Equation: Displacement = ({tunnel_length_m} / 2) + ({moth_airspeed_m_per_min} - {air_speed_m_per_min}) * ((({num_leds} - 1) * {led_delay_s}) / {seconds_per_minute})")
    print(f"Calculated Displacement = {distance_to_trigger_m}m + {moth_ground_speed_m_per_min:.3f}m/min * {pulse_travel_time_min:.3f}min")
    print(f"Total Displacement = {total_displacement_m:.2f}m")

solve_moth_displacement()