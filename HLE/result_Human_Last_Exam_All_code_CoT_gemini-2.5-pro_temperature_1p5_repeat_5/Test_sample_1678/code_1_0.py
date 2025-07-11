import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 0: Given parameters
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56  # Displayed red time (seconds)
    AR = 2      # All-red time (seconds)
    g = 30      # Effective green time (seconds)
    v_hr = 600   # Approach average flow rate (veh/hour)
    percent_green_arrivals = 0.40
    percent_red_arrivals = 0.60

    # Step 1: Calculate Cycle Parameters
    # Convert flow rates from veh/hour to veh/second
    s = s_hr / 3600
    v = v_hr / 3600
    
    # Calculate effective red time (r) and cycle length (C)
    r = R_disp + AR
    C = g + r

    print("--- Step 1: Cycle Parameters ---")
    print(f"Effective red time (r) = Displayed Red + All-Red = {R_disp}s + {AR}s = {r}s")
    print(f"Cycle length (C) = Effective Green + Effective Red = {g}s + {r}s = {C}s")
    print(f"Saturation flow rate (s) = {s_hr} veh/hr = {s:.4f} veh/s")
    print(f"Average flow rate (v) = {v_hr} veh/hr = {v:.4f} veh/s")
    print("-" * 35)

    # Step 2: Determine Arrival Characteristics
    # Total vehicles per cycle
    N = v * C
    # Vehicles arriving during red and green intervals
    N_red_arrivals = percent_red_arrivals * N
    N_green_arrivals = percent_green_arrivals * N
    
    # Arrival rates during red (lambda2) and green (lambda1)
    lambda2 = N_red_arrivals / r
    lambda1 = N_green_arrivals / g

    # The queue at the start of green is the number of vehicles that arrived during red
    Q_start_green = N_red_arrivals

    print("--- Step 2: Arrival Characteristics ---")
    print(f"Total vehicles per cycle (N) = v * C = {v:.4f} * {C} = {N:.4f} veh")
    print(f"Queue at start of green (arrivals during red) = {percent_red_arrivals:.2f} * {N:.4f} = {Q_start_green:.4f} veh")
    print(f"Arrival rate during green (λ1) = ({percent_green_arrivals:.2f} * {N:.4f}) / {g}s = {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2) = {Q_start_green:.4f} / {r}s = {lambda2:.4f} veh/s")
    print("-" * 35)
    
    # Step 3: Calculate Queue Clearance Time
    # Time to clear the queue, measured from the start of green
    t_clear_g = Q_start_green / (s - lambda1)
    
    print("--- Step 3: Queue Clearance Time ---")
    print(f"Time to clear queue (t_clear_g) = Q_start_green / (s - λ1)")
    print(f"t_clear_g = {Q_start_green:.4f} / ({s:.4f} - {lambda1:.4f}) = {t_clear_g:.2f}s")
    print(f"This is less than the effective green time of {g}s, so the queue clears.")
    print("-" * 35)

    # Step 4: Calculate Total and Average Delay
    # Total delay per cycle (W) using the area of the queuing diagram
    W = 0.5 * Q_start_green * (r + t_clear_g)
    
    # Average delay per vehicle (d)
    d = W / N
    
    print("--- Step 4: Calculate Average Delay ---")
    print(f"Total delay per cycle (W) = 0.5 * Q_start_green * (r + t_clear_g)")
    print(f"W = 0.5 * {Q_start_green:.4f} veh * ({r}s + {t_clear_g:.2f}s) = {W:.2f} veh-s")
    print("\nFinal Calculation:")
    print(f"Average Delay (d) = Total Delay (W) / Total Vehicles (N)")
    print(f"d = {W:.2f} / {N:.4f} = {d:.2f} seconds per vehicle")


calculate_deterministic_delay()
<<<23.93>>>