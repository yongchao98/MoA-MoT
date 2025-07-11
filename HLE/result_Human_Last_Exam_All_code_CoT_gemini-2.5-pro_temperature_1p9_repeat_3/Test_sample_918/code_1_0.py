import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement based on the problem description.
    """
    # Step 1: Define constants and calculate the moth's ground speed
    moth_airspeed_mpm = 5.675  # m/min
    wind_speed_mpm = 5.0      # m/min
    moth_ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm

    tunnel_length_m = 2.0
    halfway_distance_m = tunnel_length_m / 2

    # Step 2: Calculate the time for the blinking wave to travel
    num_lights = 80
    intervals = num_lights - 1
    delay_per_light_s = 0.3  # seconds
    
    # Total time for the blink wave to travel from west to east in seconds
    wave_travel_time_s = intervals * delay_per_light_s
    
    # Convert wave travel time to minutes to match the speed unit
    wave_travel_time_min = wave_travel_time_s / 60.0

    # Step 3: Calculate the additional distance traveled by the moth
    # The moth is already at the halfway point (1m) when the blinking starts.
    # We need to calculate how much farther it flies while the blinking wave travels.
    extra_distance_m = moth_ground_speed_mpm * wave_travel_time_min

    # The total displacement is the initial 1m plus the extra distance traveled.
    total_displacement_m = halfway_distance_m + extra_distance_m

    # Step 4: Print the final equation and the result
    print("The final displacement is the sum of the distance to the halfway point and the additional distance the moth flies while the LED signal travels.")
    print("\nFinal Equation:")
    print(f"Displacement = (Halfway Distance) + (Moth Ground Speed) * (Wave Travel Time)")
    
    # Showing the calculation with numbers and consistent units (meters, minutes)
    print(f"Displacement = {halfway_distance_m}m + ({moth_airspeed_mpm}m/min - {wind_speed_mpm}m/min) * (({num_lights} - 1) * {delay_per_light_s}s / 60s/min)")
    print(f"Displacement = {halfway_distance_m}m + ({moth_ground_speed_mpm}m/min) * ({wave_travel_time_min:.3f}min)")
    print(f"Displacement = {halfway_distance_m}m + {extra_distance_m:.6f}m")

    print(f"\nFinal calculated displacement: {total_displacement_m:.6f}m")
    # Rounding to match the answer choices
    print(f"Rounded to two decimal places: {round(total_displacement_m, 2)}m")


solve_moth_displacement()
<<<D>>>