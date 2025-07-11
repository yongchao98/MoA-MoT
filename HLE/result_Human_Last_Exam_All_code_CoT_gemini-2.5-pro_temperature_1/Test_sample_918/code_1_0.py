import math

def solve_moth_problem():
    """
    This function calculates the moth's displacement based on the given parameters.
    """
    # --- Givens ---
    tunnel_length = 2.0  # meters
    air_speed_east = 5.0  # m/min
    moth_air_speed_west = 5.675  # m/min
    num_leds = 80
    led_delay_sec = 0.3  # seconds

    # --- Calculations ---

    # 1. Moth's ground speed (flying west against the east wind)
    # The moth's speed relative to the ground is its flying speed minus the headwind.
    moth_ground_speed_mpm = moth_air_speed_west - air_speed_east

    # 2. Time until the light sequence starts
    # The sequence starts when the moth is halfway, i.e., has traveled 1m.
    distance_to_halfway_m = tunnel_length / 2.0
    time_to_halfway_min = distance_to_halfway_m / moth_ground_speed_mpm

    # 3. Time for the blink signal to travel the full length
    # There are 80 lights, meaning 79 intervals between them.
    num_intervals = num_leds - 1
    blink_travel_time_sec = num_intervals * led_delay_sec
    # Convert this time to minutes to be consistent.
    blink_travel_time_min = blink_travel_time_sec / 60.0

    # 4. Total time elapsed until the last LED blinks
    total_time_min = time_to_halfway_min + blink_travel_time_min

    # 5. Total displacement of the moth
    # Displacement = Speed * Time
    total_displacement_m = moth_ground_speed_mpm * total_time_min

    # --- Output ---
    print("The goal is to find the moth's total displacement from its start (the eastern end).")
    print("\nStep 1: Calculate the moth's speed relative to the ground.")
    print(f"Moth Ground Speed = Moth Air Speed - Wind Speed")
    print(f"Moth Ground Speed = {moth_air_speed_west} m/min - {air_speed_east} m/min = {moth_ground_speed_mpm} m/min")

    print("\nStep 2: Calculate the total time until the last LED blinks.")
    print("This is the time to reach the halfway point plus the time for the light signal to travel.")
    print(f"Time to halfway = (Distance to halfway) / (Ground Speed) = ({tunnel_length}/2)m / {moth_ground_speed_mpm} m/min = {time_to_halfway_min:.4f} min")
    print(f"Light signal travel time = ({num_leds} - 1) intervals * {led_delay_sec} s/interval / 60 s/min = {blink_travel_time_min:.4f} min")
    print(f"Total time = {time_to_halfway_min:.4f} min + {blink_travel_time_min:.4f} min = {total_time_min:.4f} min")

    print("\nStep 3: Calculate the final displacement.")
    print("Displacement = Ground Speed * Total Time")
    # The calculation can also be expressed as: Displacement = Distance to Halfway + (Ground Speed * Light Signal Travel Time)
    # This simplifies the equation.
    print(f"Final Equation: ({tunnel_length} / 2) + ({moth_air_speed_west} - {air_speed_east}) * (({num_leds} - 1) * {led_delay_sec} / 60) = {total_displacement_m:.2f}m")

solve_moth_problem()
<<<D>>>