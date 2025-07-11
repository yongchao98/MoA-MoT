import math

# Define the constants from the problem description
tunnel_length_m = 2.0
wind_speed_m_per_min = 5.0
moth_airspeed_m_per_min = 5.675
num_leds = 80
led_delay_s = 0.3

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the problem description.
    """
    # Step 1: Calculate the moth's ground speed
    # The moth flies against the wind, so we subtract the wind speed from the moth's airspeed.
    print("Step 1: Calculate the moth's effective speed relative to the ground.")
    moth_ground_speed_m_per_min = moth_airspeed_m_per_min - wind_speed_m_per_min
    print(f"The moth's ground speed is its airspeed ({moth_airspeed_m_per_min} m/min) minus the headwind speed ({wind_speed_m_per_min} m/min).")
    print(f"Ground Speed = {moth_airspeed_m_per_min} - {wind_speed_m_per_min} = {moth_ground_speed_m_per_min} m/min.\n")

    # Step 2: Identify the trigger point for the LED sequence
    # The sequence starts when the moth is halfway down the tunnel.
    halfway_point_m = tunnel_length_m / 2.0
    print("Step 2: Identify the trigger event for the LED blinking.")
    print(f"The LED sequence starts when the moth reaches the halfway point of the {tunnel_length_m}m tunnel, which is at {halfway_point_m}m.\n")

    # Step 3: Calculate the duration of the LED sequence
    # The time from the first LED (western end) to the last LED (eastern end).
    print("Step 3: Calculate the total time for the blink to travel from the first to the last LED.")
    num_intervals = num_leds - 1
    total_delay_s = num_intervals * led_delay_s
    print(f"There are {num_leds} LEDs, which means there are {num_leds} - 1 = {num_intervals} time intervals.")
    print(f"Total sequence time = {num_intervals} intervals * {led_delay_s} seconds/interval = {total_delay_s} seconds.\n")

    # Step 4: Calculate the distance the moth travels during the LED sequence
    print("Step 4: Calculate how far the moth travels during the LED sequence.")
    # To use the delay in seconds, we convert the moth's speed from m/min to m/s.
    moth_ground_speed_m_per_s = moth_ground_speed_m_per_min / 60.0
    distance_traveled_during_delay = moth_ground_speed_m_per_s * total_delay_s
    print(f"Distance = speed * time = ({moth_ground_speed_m_per_min}/60 m/s) * {total_delay_s} s = {distance_traveled_during_delay:.4f} m.\n")

    # Step 5: Calculate the final displacement
    print("Step 5: Calculate the final displacement from the starting point.")
    final_displacement = halfway_point_m + distance_traveled_during_delay
    print("The moth was already at the halfway point when the sequence began. We add the distance it traveled during the sequence to find its final position.")
    print(f"Final Displacement = (Position at start of sequence) + (Distance traveled during sequence)")
    print(f"Final Displacement = {halfway_point_m} m + {distance_traveled_during_delay:.4f} m = {final_displacement:.4f} m.")
    
    # Round to match the format of the answer choices
    print(f"\nRounded to two decimal places, the moth's displacement is {final_displacement:.2f}m.")

solve_moth_displacement()
<<<D>>>