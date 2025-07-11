import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection
    with non-uniform arrival rates.
    """
    # Given parameters
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56  # Displayed red time (seconds)
    Y = 3        # Displayed yellow time (seconds)
    AR = 2       # All-red time (seconds)
    g = 30       # Effective green time (seconds)
    t_L = 4      # Total lost time (seconds)
    v_hr = 600   # Average approach flow rate (veh/hour)
    percent_g_arrival = 0.4 # Percentage of traffic arriving during green
    percent_r_arrival = 0.6 # Percentage of traffic arriving during red

    # --- Step 1: Determine Signal Cycle Characteristics ---
    
    # Convert flow rates from veh/hour to veh/second
    s = s_hr / 3600
    v_avg = v_hr / 3600

    # From g = G_disp + Y + AR - t_L, we find the displayed green time (G_disp)
    G_disp = g - Y - AR + t_L
    
    # Calculate the total cycle length (C)
    C = G_disp + Y + R_disp
    
    # Calculate the effective red time (r)
    r = C - g
    
    # --- Step 2: Calculate Non-Uniform Arrival Rates ---

    # Total vehicles arriving per cycle (N)
    N = v_avg * C

    # Number of vehicles arriving during effective green (N_g) and red (N_r)
    N_g = N * percent_g_arrival
    N_r = N * percent_r_arrival
    
    # Arrival rate during effective green (lambda1)
    lambda1 = N_g / g

    # Arrival rate during effective red (lambda2)
    lambda2 = N_r / r

    # --- Step 3: Calculate Total Delay per Cycle (D_cycle) ---

    # The formula for total delay per cycle is derived from the area of the D/D/1 queueing diagram.
    # It has two parts:
    # 1. Delay of vehicles arriving during red, over the duration of red.
    delay_term1 = 0.5 * lambda2 * r**2
    # 2. Delay of all queued vehicles until the queue clears during green.
    delay_term2 = (0.5 * (lambda2 * r)**2) / (s - lambda1)
    
    D_cycle = delay_term1 + delay_term2
    
    # --- Step 4: Calculate Average Delay per Vehicle (d) ---
    
    d = D_cycle / N
    
    # --- Print Results ---
    
    print(f"The final equation for average delay is: d = (Term1 + Term2) / N")
    print("Where:")
    print(f"Term1 = Delay from red period arrivals = 0.5 * {lambda2:.4f} * {r:.2f}^2 = {delay_term1:.2f} veh-s")
    print(f"Term2 = Additional delay until queue clear = 0.5 * ({lambda2:.4f} * {r:.2f})^2 / ({s:.4f} - {lambda1:.4f}) = {delay_term2:.2f} veh-s")
    print(f"N = Total vehicles per cycle = {N:.2f} veh")
    print("\nPlugging in the numbers:")
    final_eq_str = f"d = ({delay_term1:.2f} + {delay_term2:.2f}) / {N:.2f}"
    print(final_eq_str)
    
    print(f"d = {D_cycle:.2f} / {N:.2f}")
    
    print(f"\nThe average deterministic delay per vehicle is {d:.2f} seconds.")

# Run the calculation
calculate_deterministic_delay()
print("<<<23.93>>>")