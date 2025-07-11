import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the problem's conditions.
    """
    # Given parameters
    tunnel_length_m = 2.0
    wind_speed_m_per_min = 5.0
    moth_airspeed_m_per_min = 5.675
    led_count = 80
    led_blink_delay_s = 0.3

    # Step 1: Calculate the moth's ground speed in m/s
    # The moth flies west, and the wind blows east, so their speeds subtract.
    moth_ground_speed_m_per_min = moth_airspeed_m_per_min - wind_speed_m_per_min
    # Convert speed to meters per second for consistency in calculations
    moth_ground_speed_m_per_s = moth_ground_speed_m_per_min / 60.0

    # Step 2: Calculate the total time until the eastern-most LED blinks
    
    # Part A: Calculate the time for the moth to reach the halfway point of the tunnel.
    # The first LED blink is triggered when the moth is at the 1m mark.
    halfway_distance_m = tunnel_length_m / 2.0
    time_to_halfway_s = halfway_distance_m / moth_ground_speed_m_per_s
    
    # Part B: Calculate the time for the blink signal to propagate from the first to the last LED.
    # There are (led_count - 1) intervals between the 80 LEDs.
    led_propagation_time_s = (led_count - 1) * led_blink_delay_s
    
    # The total time is the sum of the moth's travel time to halfway and the signal propagation time.
    total_time_s = time_to_halfway_s + led_propagation_time_s
    
    # Step 3: Calculate the moth's final displacement from the start
    # Displacement = Speed * Time
    final_displacement_m = moth_ground_speed_m_per_s * total_time_s
    
    # Print the final equation as requested
    print("The moth's final displacement is its ground speed multiplied by the total elapsed time.")
    print("This can be represented by the equation: Displacement = Speed * Time")
    print("\nCalculated values for the final equation:")
    # Using f-string formatting to display the numbers in the final equation.
    # f-string formatting rounds the result for display purposes, but the calculation uses the full precision value.
    print(f"{moth_ground_speed_m_per_s:.5f} m/s * {total_time_s:.3f} s = {final_displacement_m:.2f} m")

# Run the solver
solve_moth_displacement()
<<<D>>>