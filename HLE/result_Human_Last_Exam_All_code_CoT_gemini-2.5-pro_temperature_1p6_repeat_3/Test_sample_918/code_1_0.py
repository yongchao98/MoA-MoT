import math

def solve_moth_problem():
    """
    This function calculates the moth's displacement based on the given parameters.
    """
    # --- Step 1: Define constants and calculate moth's ground speed ---
    tunnel_length_m = 2.0
    wind_speed_mpm = 5.0      # meters per minute, east
    moth_airspeed_mpm = 5.675 # meters per minute, west (against the wind)
    num_leds = 80
    led_delay_s = 0.3

    # The moth's speed relative to the ground is its airspeed minus the headwind speed.
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
    # Convert speed from meters per minute to meters per second.
    moth_ground_speed_mps = moth_ground_speed_mpm / 60.0

    # --- Step 2: Calculate the blink propagation time ---
    # The blinking signal travels across (num_leds - 1) intervals.
    num_intervals = num_leds - 1
    blink_travel_time_s = num_intervals * led_delay_s

    # --- Step 3: Determine moth's position ---
    # The blinking sequence starts when the moth is at the halfway point.
    halfway_distance_m = tunnel_length_m / 2.0
    # Calculate the additional distance the moth travels during the blink propagation time.
    additional_distance_m = moth_ground_speed_mps * blink_travel_time_s

    # --- Step 4: Calculate final displacement ---
    # The final displacement is the halfway distance plus the additional distance traveled.
    final_displacement_m = halfway_distance_m + additional_distance_m

    # --- Output the results as a clear, step-by-step equation ---
    print("The moth's final displacement is its position when the blink sequence started plus the additional distance it flies while the blink signal travels.")
    print("Final Displacement = Start Position of Event + (Moth's Ground Speed * Blink Travel Time)\n")
    
    print("Breaking down the final equation:")
    print(f"Start Position of Event (Halfway Point) = {halfway_distance_m} m")
    
    # Show the calculation for the additional distance
    print(f"Additional Distance = ({moth_airspeed_mpm} m/min - {wind_speed_mpm} m/min)/60s/min * (({num_leds} - 1) * {led_delay_s} s)")
    print(f"Additional Distance = ({moth_ground_speed_mpm} m/min / 60) * ({num_intervals} * {led_delay_s}) s")
    print(f"Additional Distance = {moth_ground_speed_mps:.5f} m/s * {blink_travel_time_s} s")
    print(f"Additional Distance = {additional_distance_m:.6f} m\n")

    # Show the final calculation
    print("Final Calculation:")
    print(f"Displacement = {halfway_distance_m} m + {additional_distance_m:.6f} m")
    print(f"Displacement = {final_displacement_m:.6f} m")
    
    print(f"\nThe moth's displacement is approximately {final_displacement_m:.2f}m.")


solve_moth_problem()
<<<D>>>