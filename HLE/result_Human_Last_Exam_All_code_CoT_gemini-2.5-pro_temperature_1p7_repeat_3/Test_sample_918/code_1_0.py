import math

def solve_moth_problem():
    """
    Calculates the moth's displacement when the easternmost LED blinks.
    """
    # Define the given parameters
    tunnel_length_m = 2.0
    moth_speed_air_m_per_min = 5.675
    wind_speed_m_per_min = 5.0
    led_count = 80
    led_delay_s = 0.3

    # --- Step 1: Calculate the moth's net speed relative to the ground ---
    # The moth flies against the wind, so its net speed is its speed in air minus the wind speed.
    net_speed_m_per_min = moth_speed_air_m_per_min - wind_speed_m_per_min
    
    # Convert speed from meters per minute to meters per second for calculations involving seconds.
    secs_per_min = 60.0
    net_speed_m_per_s = net_speed_m_per_min / secs_per_min

    # --- Step 2: Determine timing and displacement components ---
    # The LED blinking starts when the moth is halfway, which is its initial position for the second phase of our calculation.
    position_at_cycle_start_m = tunnel_length_m / 2.0
    
    # Calculate the total time it takes for the blink to travel from the first (western) LED to the last (eastern) LED.
    # There are (led_count - 1) intervals between the 80 LEDs.
    led_propagation_time_s = (led_count - 1) * led_delay_s
    
    # --- Step 3: Calculate the final displacement ---
    # During the LED propagation time, the moth continues to move.
    # Calculate the additional distance traveled during this time.
    additional_distance_m = net_speed_m_per_s * led_propagation_time_s
    
    # The final displacement is its position when the cycle started plus the additional distance it traveled.
    final_displacement_m = position_at_cycle_start_m + additional_distance_m

    # --- Output the final equation and result ---
    print("The moth's final displacement is its position when the LED blinking sequence starts, plus the additional distance it travels while the blinking signal propagates from the west to the east end.")
    print("\nFinal Displacement = Position_at_Cycle_Start + (Net_Speed * LED_Propagation_Time)")
    
    print("\nBreaking down the equation:")
    print(f"Position_at_Cycle_Start = {position_at_cycle_start_m} m")
    
    # Note: For clarity, Net_Speed is shown in m/s, calculated from m/min.
    print(f"Net_Speed = ({moth_speed_air_m_per_min} - {wind_speed_m_per_min}) / {secs_per_min} = {net_speed_m_per_s:.5f} m/s")
    print(f"LED_Propagation_Time = ({led_count} - 1) * {led_delay_s} = {led_propagation_time_s:.1f} s")
    
    print("\nFinal Equation with values:")
    print(f"Final Displacement = {position_at_cycle_start_m} m + ({net_speed_m_per_s:.5f} m/s * {led_propagation_time_s:.1f} s)")
    print(f"Final Displacement = {position_at_cycle_start_m} m + {additional_distance_m:.5f} m")
    print(f"Final Displacement = {final_displacement_m:.5f} m")
    print(f"\nThe moth's displacement, rounded to two decimal places, is {final_displacement_m:.2f} m.")


solve_moth_problem()
<<<D>>>