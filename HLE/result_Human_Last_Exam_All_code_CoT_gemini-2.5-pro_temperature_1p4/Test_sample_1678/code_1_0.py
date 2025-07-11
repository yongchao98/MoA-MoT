import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 1: Define given parameters
    s_vph = 2160.0  # Saturation flow rate (veh/hour)
    R = 56.0        # Displayed red time (s)
    Y = 3.0         # Displayed yellow time (s)
    AR = 2.0        # All-red time (s)
    g = 30.0        # Effective green time (s)
    v_vph = 600.0   # Approach average flow rate (veh/hour)
    percent_red_arrivals = 0.60   # 60% of traffic arrives during red
    percent_green_arrivals = 0.40 # 40% of traffic arrives during green

    # Step 2: Convert units and calculate cycle parameters
    s = s_vph / 3600.0  # Saturation flow rate (veh/s)
    v = v_vph / 3600.0  # Average arrival rate (veh/s)

    # Effective red time (r) is the period when vehicles are not served.
    r = R + AR
    # Cycle length (C) is the sum of effective green and effective red times.
    C = g + r

    print("--- Step 1 & 2: Initial Parameters and Cycle Characteristics ---")
    print(f"Saturation flow rate (s): {s_vph} veh/h = {s:.4f} veh/s")
    print(f"Average arrival rate (v): {v_vph} veh/h = {v:.4f} veh/s")
    print(f"Effective red time (r) = Displayed Red + All-Red = {R} + {AR} = {r} s")
    print(f"Cycle length (C) = Effective Green + Effective Red = {g} + {r} = {C} s\n")

    # Step 3: Determine arrival profiles
    total_arrivals_per_cycle = v * C
    arrivals_during_red = total_arrivals_per_cycle * percent_red_arrivals
    arrivals_during_green = total_arrivals_per_cycle * percent_green_arrivals

    # Arrival rate during the effective red interval
    lambda2 = arrivals_during_red / r
    # Arrival rate during the effective green interval
    lambda1 = arrivals_during_green / g

    print("--- Step 3: Arrival Profile Calculation ---")
    print(f"Total arrivals per cycle (N) = v * C = {v:.4f} * {C} = {total_arrivals_per_cycle:.4f} vehicles")
    print(f"Arrival rate during red (λ2) = ({total_arrivals_per_cycle:.4f} * {percent_red_arrivals}) / {r} = {lambda2:.4f} veh/s")
    print(f"Arrival rate during green (λ1) = ({total_arrivals_per_cycle:.4f} * {percent_green_arrivals}) / {g} = {lambda1:.4f} veh/s\n")

    # Step 4: Calculate total delay per cycle (D_total)
    # The maximum queue at the end of the red interval is the number of vehicles that arrived during red.
    Q_max = arrivals_during_red

    # The time it takes for this queue to clear from the start of green.
    # The queue dissipates at a rate of (s - λ1).
    time_to_clear_queue = Q_max / (s - lambda1)

    # The total delay is the area under the queue length curve over one cycle.
    # This can be calculated as the sum of delay during red and delay during green clearance.
    # Delay during red = Area of a triangle with base r and height Q_max
    delay_during_red = 0.5 * Q_max * r
    # Delay during green clearance = Area of a triangle with base t_c and height Q_max
    delay_during_green_clearance = 0.5 * Q_max * time_to_clear_queue

    total_delay_per_cycle = delay_during_red + delay_during_green_clearance
    
    print("--- Step 4: Total Delay Calculation ---")
    print(f"Maximum queue at start of green (Q_max) = {Q_max:.4f} vehicles")
    print(f"Time to clear queue (t_c) = Q_max / (s - λ1) = {Q_max:.4f} / ({s:.4f} - {lambda1:.4f}) = {time_to_clear_queue:.4f} s")
    print(f"Total Delay per Cycle (D_total) = (0.5 * Q_max * r) + (0.5 * Q_max * t_c)")
    print(f"D_total = (0.5 * {Q_max:.4f} * {r}) + (0.5 * {Q_max:.4f} * {time_to_clear_queue:.4f})")
    print(f"D_total = {delay_during_red:.4f} + {delay_during_green_clearance:.4f} = {total_delay_per_cycle:.4f} veh-s\n")
    
    # Step 5: Calculate average delay per vehicle
    average_delay = total_delay_per_cycle / total_arrivals_per_cycle

    print("--- Step 5: Average Delay Per Vehicle ---")
    print(f"Average Delay (d) = Total Delay / Total Arrivals")
    print(f"d = {total_delay_per_cycle:.4f} / {total_arrivals_per_cycle:.4f}")
    print(f"The average deterministic delay per vehicle is {average_delay:.2f} seconds.")

    return average_delay

# Run the calculation and print the final answer in the required format
final_answer = calculate_deterministic_delay()
print(f"<<<{final_answer:.2f}>>>")
