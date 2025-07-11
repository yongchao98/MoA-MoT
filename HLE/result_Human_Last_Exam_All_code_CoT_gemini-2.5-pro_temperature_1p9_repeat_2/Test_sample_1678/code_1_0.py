import math

def solve_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a pre-timed signalized intersection.
    """
    # 1. Define input variables
    s_h = 2160  # veh/hour
    R_disp = 56 # s
    Y = 3       # s
    AR = 2      # s
    g = 30      # s
    t_L = 4     # s
    v_h = 600   # veh/hour
    p_green_arrival = 0.40 # percentage of traffic arriving during green
    p_red_arrival = 0.60   # percentage of traffic arriving during red

    # --- Calculations ---

    # Convert flow rates from vehicles per hour to vehicles per second
    s = s_h / 3600
    v = v_h / 3600

    # Step 1: Calculate Cycle Length (C)
    # Effective green g = G_disp + Y + AR - t_L => G_disp = g - Y - AR + t_L
    G_disp = g - Y - AR + t_L
    # Cycle length C = G_disp + Y + R_disp
    C = G_disp + Y + R_disp

    # Step 2: Calculate Effective Red Time (r)
    # Effective red r is the part of the cycle that is not effective green
    r = C - g

    # Step 3: Calculate arrivals per cycle
    # Total vehicles arriving in one cycle
    N = v * C
    # Vehicles arriving during the effective red interval
    N_r = p_red_arrival * N
    # Vehicles arriving during the effective green interval
    N_g = p_green_arrival * N

    # Step 4: Calculate time to clear the queue that built during red (t_clear)
    # The arrival rate during the effective green interval (λ1)
    lambda_1 = N_g / g
    
    # The time to clear is the queue at the end of red (N_r) divided by the
    # rate at which the queue is processed (s - λ1)
    if s <= lambda_1:
        print("Error: Saturation flow rate is not greater than arrival rate during green.")
        print("The queue will never clear. The intersection is oversaturated.")
        return

    t_clear = N_r / (s - lambda_1)
    
    # --- Output ---

    print("This program calculates the average deterministic delay per vehicle (d) using D/D/1 queuing principles with non-uniform arrival rates.\n")

    print("--- Intermediate Calculations ---\n")

    print("Step 1: Calculate Cycle Length (C)")
    print(f"First, we find the displayed green time (G_disp):")
    print(f"G_disp = g - Y - AR + t_L = {g} - {Y} - {AR} + {t_L} = {G_disp:.2f} s")
    print("Now, we calculate the cycle length (C):")
    print(f"C = G_disp + Y + R_disp = {G_disp:.2f} + {Y} + {R_disp} = {C:.2f} s\n")

    print("Step 2: Calculate Effective Red Time (r)")
    print(f"r = C - g = {C:.2f} - {g} = {r:.2f} s\n")

    print("Step 3: Calculate Vehicle Arrivals per Cycle")
    print(f"Total vehicles per cycle (N) = Average flow rate (v) * C")
    print(f"N = {v:.4f} veh/s * {C:.2f} s = {N:.4f} vehicles")
    print(f"Vehicles arriving during red (N_r) = {p_red_arrival:.2f} * {N:.4f} = {N_r:.4f} vehicles")
    print(f"Vehicles arriving during green (N_g) = {p_green_arrival:.2f} * {N:.4f} = {N_g:.4f} vehicles\n")

    print("Step 4: Calculate Time to Clear the Red-Time Queue (t_clear)")
    print(f"Arrival rate during green (λ1) = N_g / g = {N_g:.4f} / {g} = {lambda_1:.4f} veh/s")
    print(f"Saturation flow rate (s) = {s_h} veh/hr / 3600 = {s:.4f} veh/s")
    print(f"t_clear = N_r / (s - λ1) = {N_r:.4f} / ({s:.4f} - {lambda_1:.4f}) = {t_clear:.4f} s")

    if t_clear > g:
        print("Warning: The queue clearance time is longer than the effective green time.")
        print("This indicates an oversaturated condition where not all queue clears each cycle.\n")
    else:
        print("This is less than the effective green time, so the queue clears each cycle.\n")

    # Step 5: Calculate Total Delay per Cycle (TD_cycle)
    # This is the area of the queuing diagram, composed of delay during red and delay during green while the queue clears.
    # TD_cycle = (0.5 * N_r * r) + (0.5 * N_r * t_clear)
    TD_cycle = 0.5 * N_r * (r + t_clear)
    
    print("Step 5: Calculate Total Delay per Cycle (TD_cycle)")
    print("TD_cycle = 0.5 * N_r * (r + t_clear)")
    print(f"TD_cycle = 0.5 * {N_r:.4f} * ({r:.2f} + {t_clear:.4f}) = {TD_cycle:.4f} veh-seconds\n")
    
    # Step 6: Calculate Average Delay per Vehicle (d)
    d = TD_cycle / N

    print("--- Final Result ---\n")
    print("Step 6: Calculate Average Delay per Vehicle (d)")
    print("d = Total Delay per Cycle (TD_cycle) / Total Vehicles per Cycle (N)")
    print(f"d = {TD_cycle:.4f} / {N:.4f} = {d:.2f} seconds\n")

if __name__ == '__main__':
    solve_traffic_delay()

<<<23.93>>>