import math

def solve_moth_problem():
    """
    This function calculates the moth's displacement based on the provided scenario.
    """
    # Step 1: Define variables and convert units
    tunnel_length_m = 2.0
    wind_speed_m_per_min = 5.0
    moth_airspeed_m_per_min = 5.675
    num_leds = 80
    led_delay_s = 0.3
    
    # Convert speeds from m/min to m/s
    moth_airspeed_m_per_s = moth_airspeed_m_per_min / 60
    wind_speed_m_per_s = wind_speed_m_per_min / 60
    
    # Step 2: Calculate the moth's effective ground speed
    # The moth flies west, against the eastward wind.
    moth_ground_speed_m_per_s = moth_airspeed_m_per_s - wind_speed_m_per_s
    
    # Step 3: Define constants for the displacement calculation
    # The LED blinking starts when the moth has traveled 1m from the eastern end.
    distance_to_halfway_m = tunnel_length_m / 2.0
    
    # The blinking propagates from the 1st to the 80th LED, which is 79 steps.
    num_led_gaps = num_leds - 1
    
    # The displacement is the total distance the moth has travelled.
    # This can be calculated directly as:
    # Displacement = (Distance to halfway) + (Distance travelled while LEDs blink)
    # Distance travelled while LEDs blink = moth_ground_speed * time_for_leds_to_blink
    # time_for_leds_to_blink = num_led_gaps * led_delay
    
    displacement_m = distance_to_halfway_m + (moth_ground_speed_m_per_s * num_led_gaps * led_delay_s)
    
    # Step 4: Output the calculation and the result
    print("The moth's displacement is the distance it travels to the halfway point plus the distance it travels while the LED signal propagates.")
    print("\nFinal Equation:")
    print(f"Displacement = Distance to Halfway + (Moth Ground Speed * (Number of LEDs - 1) * LED Delay)")
    print("\nPlugging in the numbers:")
    # We round the ground speed for display purposes but use the full precision in calculation.
    print(f"{displacement_m:.2f}m = {distance_to_halfway_m}m + ({moth_ground_speed_m_per_s:.5f}m/s * {num_led_gaps} * {led_delay_s}s)")
    
solve_moth_problem()
<<<D>>>