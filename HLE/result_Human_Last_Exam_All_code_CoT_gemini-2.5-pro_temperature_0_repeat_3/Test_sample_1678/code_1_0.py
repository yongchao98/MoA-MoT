import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # --- Given Information ---
    s_hr = 2160      # Saturation flow rate (veh/hour)
    R = 56           # Displayed red time (seconds)
    Y = 3            # Displayed yellow time (seconds)
    AR = 2           # All-red time (seconds)
    g = 30           # Effective green time (seconds)
    tL = 4           # Total lost time (seconds)
    v_hr = 600       # Average approach flow rate (veh/hour)
    pct_green_arrival = 0.40 # Percentage of traffic arriving during green
    pct_red_arrival = 0.60   # Percentage of traffic arriving during red

    # --- Step 1: Calculate Cycle Time and Effective Red ---
    print("--- Step 1: Calculate Cycle Characteristics ---")
    # The effective green time is related to the displayed times and lost time by: g = G + Y + AR - tL
    # We can find the displayed green time G from this formula.
    G = g - Y - AR + tL
    print(f"1a. Calculate Displayed Green (G):")
    print(f"    G = g - Y - AR + tL = {g}s - {Y}s - {AR}s + {tL}s = {G:.2f} s")

    # The cycle length C is the sum of the displayed signal times.
    C = G + Y + R
    print(f"1b. Calculate Cycle Length (C):")
    print(f"    C = G + Y + R = {G:.2f}s + {Y}s + {R}s = {C:.2f} s")

    # The effective red time r is the portion of the cycle that is not effectively green.
    r = C - g
    print(f"1c. Calculate Effective Red (r):")
    print(f"    r = C - g = {C:.2f}s - {g}s = {r:.2f} s\n")

    # --- Step 2: Calculate Flow Rates and Arrivals per Cycle ---
    print("--- Step 2: Calculate Flow and Arrival Rates ---")
    # Convert hourly rates to per-second rates
    s = s_hr / 3600.0
    v = v_hr / 3600.0
    print(f"2a. Convert hourly rates to veh/sec:")
    print(f"    Saturation flow rate (s) = {s_hr} veh/hr / 3600 s/hr = {s:.4f} veh/s")
    print(f"    Average arrival rate (v) = {v_hr} veh/hr / 3600 s/hr = {v:.4f} veh/s")

    # Calculate arrivals per cycle
    N = v * C
    N_r = N * pct_red_arrival
    N_g = N * pct_green_arrival
    print(f"2b. Calculate arrivals per cycle:")
    print(f"    Total arrivals per cycle (N) = {v:.4f} veh/s * {C:.2f}s = {N:.4f} veh")
    print(f"    Arrivals during red (N_r) = {N:.4f} veh * {pct_red_arrival:.2f} = {N_r:.4f} veh")
    print(f"    Arrivals during green (N_g) = {N:.4f} veh * {pct_green_arrival:.2f} = {N_g:.4f} veh")

    # Calculate arrival rates during red (lambda2) and green (lambda1)
    lambda2 = N_r / r
    lambda1 = N_g / g
    print(f"2c. Calculate specific arrival rates (λ):")
    print(f"    Arrival rate during red (λ2) = {N_r:.4f} veh / {r:.2f}s = {lambda2:.4f} veh/s")
    print(f"    Arrival rate during green (λ1) = {N_g:.4f} veh / {g:.2f}s = {lambda1:.4f} veh/s\n")

    # --- Step 3: Analyze the Queue ---
    print("--- Step 3: Analyze the Queue ---")
    # The queue at the start of green is the number of vehicles that arrived during red.
    queue_at_green_start = N_r
    # The time it takes to clear this queue depends on the dissipation rate (s - lambda1).
    time_to_clear = queue_at_green_start / (s - lambda1)
    print(f"3a. Calculate time to clear the queue after green starts:")
    print(f"    Time to clear = N_r / (s - λ1) = {queue_at_green_start:.4f} / ({s:.4f} - {lambda1:.4f}) = {time_to_clear:.4f} s")
    
    # The time in the cycle when the queue clears.
    t_c = r + time_to_clear
    print(f"3b. Calculate time in cycle when queue clears (t_c):")
    print(f"    t_c = r + time_to_clear = {r:.2f}s + {time_to_clear:.4f}s = {t_c:.4f} s\n")

    # --- Step 4: Calculate Total and Average Delay ---
    print("--- Step 4: Calculate Total and Average Delay ---")
    # Total delay is the area between the arrival and departure curves.
    # Area = Area under Arrival Curve - Area under Departure Curve
    # Area under arrival curve from 0 to t_c
    area_under_A = (0.5 * lambda2 * r**2) + (queue_at_green_start * time_to_clear + 0.5 * lambda1 * time_to_clear**2)
    # Area under departure curve from r to t_c
    area_under_D = 0.5 * s * time_to_clear**2
    
    total_delay = area_under_A - area_under_D
    print(f"4a. Calculate total delay per cycle:")
    print(f"    Total Delay = Area_under_Arrival - Area_under_Departure")
    print(f"    Total Delay = (0.5*{lambda2:.4f}*{r:.2f}² + ({N_r:.4f}*{time_to_clear:.4f} + 0.5*{lambda1:.4f}*{time_to_clear:.4f}²)) - (0.5*{s:.4f}*{time_to_clear:.4f}²)")
    print(f"    Total Delay = {area_under_A:.4f} veh-s - {area_under_D:.4f} veh-s = {total_delay:.4f} veh-s")

    # Average delay per vehicle
    avg_delay = total_delay / N
    print(f"4b. Calculate average delay per vehicle:")
    print(f"    Average Delay = Total Delay / N = {total_delay:.4f} veh-s / {N:.4f} veh = {avg_delay:.4f} s\n")

    print("--- Final Answer ---")
    print(f"The average deterministic delay per vehicle is {avg_delay:.2f} seconds.")
    return avg_delay

# Execute the function and store the final answer
final_answer = calculate_traffic_delay()
# The final answer is printed by the function, but we also return it for potential further use.
# The required format is <<<answer>>> at the end.
print(f"<<<{final_answer:.2f}>>>")
