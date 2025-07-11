import math

def calculate_average_delay():
    """
    Calculates the average deterministic delay per vehicle for a signalized intersection
    with non-uniform arrival rates.
    """
    # Step 1: Define given parameters
    s_hr = 2160.0  # Saturation flow rate (veh/hour)
    R_disp = 56.0 # Displayed red time (seconds)
    Y = 3.0       # Displayed yellow time (seconds)
    AR = 2.0      # All-red time (seconds)
    g = 30.0      # Effective green time (seconds)
    tL = 4.0      # Total lost time (seconds)
    v_hr = 600.0  # Average approach flow rate (veh/hour)
    
    # Arrival distribution
    green_arrival_pct = 0.40
    red_arrival_pct = 0.60

    # Step 2: Determine Signal Timing
    # Using the standard formula: g = G_disp + Y + AR - tL
    # We can find the displayed green time, G_disp.
    G_disp = g - Y - AR + tL
    
    # The cycle length C is the sum of the displayed intervals for the approach.
    C = G_disp + Y + R_disp
    
    # The effective red time r is the portion of the cycle that is not effective green.
    r = C - g

    # Step 3: Convert Units
    # Convert flow rates from veh/hour to veh/sec for consistency
    s = s_hr / 3600.0
    v = v_hr / 3600.0

    # Step 4: Calculate Arrival Rates (λ1, λ2)
    # Total number of vehicles arriving per cycle
    N_total_per_cycle = v * C
    
    # Number of vehicles arriving during effective green and red intervals
    N_green_arrivals = N_total_per_cycle * green_arrival_pct
    N_red_arrivals = N_total_per_cycle * red_arrival_pct
    
    # Arrival rate during effective green (λ1)
    lambda1 = N_green_arrivals / g
    
    # Arrival rate during effective red (λ2)
    lambda2 = N_red_arrivals / r

    # Step 5: Calculate Total Delay per Cycle (TD)
    # The queueing diagram method for D/D/1 with non-uniform arrivals is used.
    # Total delay is the area between the arrival and departure curves.
    
    # Delay during the effective red period (Area of first triangle in queue diagram)
    # This is the integral of the arrival function A(t) = λ2*t from 0 to r.
    delay_during_red = 0.5 * lambda2 * r**2
    
    # To calculate delay during green, first find the time to clear the queue (t_clear_green)
    # This is the time after the green starts when the departure curve catches up to the arrival curve.
    # Queue at end of red = lambda2 * r
    # This queue clears when (s - lambda1) * t_clear_green = lambda2 * r
    # Note: s > lambda1 must be true for the queue to clear.
    queue_at_end_of_red = lambda2 * r
    t_clear_green = queue_at_end_of_red / (s - lambda1)
    
    # Delay during the green period until the queue clears
    # This is the area of the second triangle in the queue diagram.
    # Base = queue_at_end_of_red, Height = t_clear_green
    delay_during_green_clearance = 0.5 * queue_at_end_of_red * t_clear_green
    
    # Total Delay per cycle
    total_delay_per_cycle = delay_during_red + delay_during_green_clearance

    # Step 6: Calculate Average Delay per Vehicle
    avg_delay_per_vehicle = total_delay_per_cycle / N_total_per_cycle
    
    # Step 7: Print the results and the final equation breakdown
    print("--- Intermediate Calculations ---")
    print(f"Displayed Green Time (G_disp): {G_disp:.2f} s")
    print(f"Cycle Length (C): {C:.2f} s")
    print(f"Effective Red Time (r): {r:.2f} s")
    print(f"Saturation Flow Rate (s): {s:.4f} veh/s")
    print(f"Average Arrival Rate (v): {v:.4f} veh/s")
    print(f"Arrival Rate during Green (λ1): {lambda1:.4f} veh/s")
    print(f"Arrival Rate during Red (λ2): {lambda2:.4f} veh/s")
    print(f"Total Vehicles per Cycle (N): {N_total_per_cycle:.2f} vehicles")
    print(f"Total Delay per Cycle: {total_delay_per_cycle:.2f} veh-s")
    
    print("\n--- Final Calculation for Average Delay ---")
    equation_str = (
        f"Average Delay = Total Delay / Total Vehicles\n"
        f"Average Delay = (Delay during Red + Delay during Green Queue Clearance) / Total Vehicles\n"
        f"Average Delay = ({delay_during_red:.2f} + {delay_during_green_clearance:.2f}) veh-s / {N_total_per_cycle:.2f} vehicles\n"
        f"Average Delay = {total_delay_per_cycle:.2f} veh-s / {N_total_per_cycle:.2f} vehicles = {avg_delay_per_vehicle:.2f} seconds"
    )
    print(equation_str)
    
    return avg_delay_per_vehicle

# Run the calculation and store the final answer
final_answer = calculate_average_delay()

# Print the final answer in the required format
print(f"\n<<<_final_answer_{final_answer:.2f}>>>")