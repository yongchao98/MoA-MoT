def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    The model is D/D/1 with different arrival rates for red and green intervals.
    """
    # Step 1: Define given parameters
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56   # Displayed red time (seconds)
    Y = 3         # Displayed yellow time (seconds)
    AR = 2        # All-red time (seconds)
    g = 30        # Effective green time (seconds)
    tL = 4        # Total lost time (seconds)
    v_avg_hr = 600 # Approach average flow rate (veh/hour)
    percent_arr_green = 0.40 # Percentage of traffic arriving during green
    percent_arr_red = 0.60   # Percentage of traffic arriving during red

    print("--- Step 1: Calculating Cycle Timing ---")
    # Calculate Displayed Green (G_disp) from the effective green time formula
    # g = G_disp + Y + AR - tL
    G_disp = g - Y - AR + tL
    print(f"Calculated Displayed Green Time (G_disp): {G_disp:.2f} s")

    # Calculate Cycle Length (C)
    # C = G_disp + Y + R_disp
    C = G_disp + Y + R_disp
    print(f"Calculated Cycle Length (C): {C:.2f} s")

    # Calculate Effective Red (r)
    # r = C - g
    r = C - g
    print(f"Calculated Effective Red Time (r): {r:.2f} s\n")

    print("--- Step 2: Calculating Rates in veh/sec ---")
    # Convert hourly rates to per-second rates
    s = s_hr / 3600
    v_avg = v_avg_hr / 3600
    print(f"Saturation Flow Rate (s): {s:.4f} veh/s")
    print(f"Average Arrival Rate (v_avg): {v_avg:.4f} veh/s\n")
    
    print("--- Step 3: Determining Arrival Characteristics ---")
    # Calculate total vehicles per cycle (N)
    N = v_avg * C
    print(f"Total Vehicles per Cycle (N): {N:.4f} veh")

    # Calculate vehicles arriving during red and green
    N_r = N * percent_arr_red
    N_g = N * percent_arr_green
    print(f"Vehicles arriving during red (N_r): {N_r:.4f} veh")
    print(f"Vehicles arriving during green (N_g): {N_g:.4f} veh")
    
    # Calculate arrival rates for red (lambda2) and green (lambda1)
    lambda2 = N_r / r if r > 0 else 0
    lambda1 = N_g / g if g > 0 else 0
    print(f"Arrival Rate during red (λ2): {lambda2:.4f} veh/s")
    print(f"Arrival Rate during green (λ1): {lambda1:.4f} veh/s\n")
    
    print("--- Step 4: Analyzing Queue Dynamics ---")
    # Maximum queue is the number of vehicles that arrive during red
    Q_max = N_r
    print(f"Maximum Queue Length (Q_max): {Q_max:.4f} veh")
    
    # Rate at which the queue dissipates
    dissipation_rate = s - lambda1
    print(f"Queue Dissipation Rate (s - λ1): {dissipation_rate:.4f} veh/s")
    
    # Time to clear the queue
    t_clear = Q_max / dissipation_rate if dissipation_rate > 0 else float('inf')
    print(f"Time to Clear Queue (t_clear): {t_clear:.4f} s\n")
    
    print("--- Step 5: Calculating Total and Average Delay ---")
    # Calculate total delay per cycle (W_c)
    # This is the area of the queueing diagram
    W_c = (0.5 * Q_max * r) + (0.5 * (Q_max**2) / dissipation_rate)
    print(f"Total Delay per Cycle (W_c): {W_c:.4f} veh-s")
    
    # Calculate average delay per vehicle (d_avg)
    d_avg = W_c / N if N > 0 else 0
    print(f"\n--- Final Calculation for Average Delay per Vehicle ---")
    print(f"Average Delay = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"              = {W_c:.2f} veh-s / {N:.2f} veh")
    print(f"              = {d_avg:.2f} s")

    return d_avg

# Execute the calculation and store the final answer
average_delay = calculate_deterministic_delay()

# Output the final answer in the required format
print(f"\nFinal Answer: The average deterministic delay per vehicle is {average_delay:.2f} seconds.")
print(f"<<<{average_delay:.2f}>>>")