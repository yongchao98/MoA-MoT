import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    based on the provided parameters.
    """
    # Step 1: Define initial parameters from the problem statement
    s_h = 2160  # Saturation flow rate (veh/hour)
    R = 56      # Displayed red time (seconds)
    Y = 3       # Displayed yellow time (seconds)
    AR = 2      # All-red time (seconds)
    g = 30      # Effective green time (seconds)
    t_L = 4     # Total lost time (seconds)
    v_h = 600   # Approach average flow rate (veh/hour)
    p_g = 0.40  # Percentage of traffic arriving during effective green
    p_r = 0.60  # Percentage of traffic arriving during effective red

    # Step 2: Determine signal timing parameters
    # Using the standard formula: g = G + Y + AR - t_L
    G = g - Y - AR + t_L
    # The cycle length C is the sum of the displayed intervals for the approach
    C = G + Y + R
    # The effective red time r is the part of the cycle that is not effective green
    r = C - g

    # Step 3: Calculate arrival and departure rates in veh/second
    s = s_h / 3600  # Saturation flow rate in veh/s
    v = v_h / 3600  # Average arrival rate in veh/s

    # Calculate total vehicles per cycle
    V_c = v * C
    
    # Calculate vehicles arriving during green and red intervals
    N_g = V_c * p_g
    N_r = V_c * p_r # This is also the maximum queue size at the end of red

    # Calculate arrival rates during green (lambda1) and red (lambda2)
    lambda1 = N_g / g
    lambda2 = N_r / r

    # Step 4: Analyze the queue to find total delay per cycle
    # Time required to clear the queue (N_r) that accumulated during the red interval.
    # The queue clears at a rate of (s - lambda1)
    if s <= lambda1:
        print("Error: Saturation flow rate is not greater than arrival rate during green.")
        print("Queue will not clear under these D/D/1 assumptions.")
        return
        
    t_clear = N_r / (s - lambda1)

    # Calculate total delay per cycle (D_total) as the area under the queueing diagram.
    # Delay accumulated during effective red time (Area 1)
    delay_during_red = 0.5 * lambda2 * r**2
    
    # Delay accumulated during effective green time until queue clears (Area 2)
    # This is the area of a triangle with base t_clear and height N_r
    delay_during_green = 0.5 * N_r * t_clear

    D_total = delay_during_red + delay_during_green

    # Step 5: Calculate the average delay per vehicle
    d = D_total / V_c
    
    # Step 6: Print the results
    print("--- Traffic Analysis Results ---")
    print(f"Cycle Length (C): {C:.2f} s")
    print(f"Effective Red Time (r): {r:.2f} s")
    print(f"Total Vehicles per Cycle (V_c): {V_c:.2f} veh")
    print(f"Arrival Rate during Red (λ2): {lambda2:.4f} veh/s")
    print(f"Arrival Rate during Green (λ1): {lambda1:.4f} veh/s")
    print(f"Queue at end of Red (N_r): {N_r:.2f} veh")
    print(f"Time to Clear Queue (t_clear): {t_clear:.2f} s")
    print(f"Total Delay per Cycle (D_total): {D_total:.2f} veh-s\n")
    
    print("--- Final Calculation ---")
    print(f"The average deterministic delay per vehicle (d) is calculated as:")
    print(f"d = Total Delay / Total Vehicles")
    print(f"d = {D_total:.2f} / {V_c:.2f}")
    print(f"d = {d:.2f} seconds\n")

# Run the calculation
calculate_deterministic_delay()
print("<<<23.95>>>")