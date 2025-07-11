import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with non-uniform arrival rates.
    """
    # --- Step 0: Define Given Parameters ---
    s_h = 2160      # Saturation flow rate (veh/hour)
    v_avg_h = 600   # Average approach flow rate (veh/hour)
    R_d = 56        # Displayed red time (s)
    Y = 3           # Displayed yellow time (s)
    AR = 2          # All-red time (s)
    g = 30          # Effective green time (s)
    t_L = 4         # Total lost time (s)
    perc_arr_green = 0.40  # Percentage of traffic arriving during effective green
    perc_arr_red = 0.60    # Percentage of traffic arriving during effective red

    # --- Step 1: Calculate Cycle Parameters ---
    # The relationship for effective green time is g = G_d + Y + AR - t_L
    # We can solve for the displayed green time, G_d.
    G_d = g - Y - AR + t_L
    # The cycle length C is the sum of all displayed time intervals.
    C = G_d + Y + R_d
    # The effective red time r is the part of the cycle that is not effective green.
    r = C - g

    # --- Step 2: Calculate Arrival Rates ---
    # Convert hourly rates to per-second rates
    s = s_h / 3600
    v_avg = v_avg_h / 3600
    # Total vehicles arriving per cycle
    N = v_avg * C
    # Number of vehicles arriving during effective green (N_g) and red (N_r)
    N_g = N * perc_arr_green
    N_r = N * perc_arr_red
    # Arrival rate during effective green (lambda1) and red (lambda2)
    lambda1 = N_g / g
    lambda2 = N_r / r

    # --- Step 3: Analyze the Queue ---
    # The queue at the end of the effective red period (Q_r) is the number of vehicles
    # that arrived during the red period.
    Q_r = N_r
    # Calculate the time it takes for this queue to clear once the green starts.
    # The queue shrinks at a rate of (s - lambda1).
    t_clear_g = Q_r / (s - lambda1)

    # --- Step 4: Calculate Total Delay per Cycle ---
    # The total delay is the area of the queue-vs-time graph.
    # Delay during red is the area of a triangle with base r and height Q_r.
    delay_red = 0.5 * Q_r * r
    # Delay during green (until queue clears) is another triangle with base t_clear_g and height Q_r.
    delay_green = 0.5 * Q_r * t_clear_g
    # Total delay per cycle is the sum of these two areas.
    D_total = delay_red + delay_green

    # --- Step 5: Calculate Average Delay per Vehicle ---
    # Average delay (d) is total delay per cycle divided by total vehicles per cycle.
    d = D_total / N

    # --- Final Output ---
    print("The calculation for the average delay per vehicle is broken down as follows:")
    print(f"Average Delay = Total Delay per Cycle / Total Vehicles per Cycle")
    # Outputting the numbers in the final equation as requested
    print(f"              = {D_total:.2f} veh-sec / {N:.2f} veh")
    print(f"              = {d:.2f} seconds")

# Execute the function to get the answer
calculate_traffic_delay()