import numpy as np

def run_simulation(E_total, S_initial, Km, kcat, total_time, time_steps):
    """
    Simulates an enzyme reaction over time.
    Returns a list of (time, product) tuples.
    """
    Vmax = kcat * E_total
    S = S_initial
    P = 0.0
    dt = total_time / float(time_steps)
    
    results = []
    time_points = np.linspace(0, total_time, time_steps + 1)
    
    # Store initial state
    results.append((0.0, 0.0))
    
    # Use a copy of time_points for iteration to avoid issues with dt
    for t in time_points[1:]:
        # Calculate reaction velocity for the current substrate concentration
        velocity = (Vmax * S) / (Km + S)
        
        # Calculate change in product and substrate
        delta_P = velocity * dt
        
        # Update concentrations for the next time step
        S -= delta_P
        P += delta_P
        
        # Ensure concentrations don't go below zero
        S = max(0, S)
        P = min(P, S_initial)
        
        results.append((t, P))
        
    return results

def main():
    """Main function to run simulations and print results."""
    # --- Enzyme and reaction parameters ---
    Km = 20.0  # Michaelis constant (uM)
    kcat = 100.0 # Turnover number (1/s)
    S_initial = 100.0 # Initial Substrate (uM)
    
    # --- Experimental conditions ---
    E_high = 1.0  # High enzyme concentration (problem case)
    E_low = 0.1   # Low enzyme concentration (solution case)
    
    total_time = 10.0 # seconds
    time_steps = 10 # number of steps in the simulation
    
    print("This script demonstrates why decreasing enzyme concentration helps reveal the linear phase of a reaction.\n")

    # --- 1. Simulate the problematic high enzyme concentration ---
    print("--- SCENARIO 1: High Enzyme Concentration (E = {:.1f} uM) ---".format(E_high))
    print("Observation: The product forms rapidly, and the rate quickly slows down (non-linear).")
    high_e_results = run_simulation(E_high, S_initial, Km, kcat, total_time, time_steps)
    
    # Print a few time points to show the non-linear curve
    # Rate = delta_P / delta_t
    print("Time (s) | Product (uM) | Rate (uM/s) in previous interval")
    print("-------- | -------------- | -----------------------------------")
    print("{:6.1f} | {:14.2f} | {}".format(high_e_results[0][0], high_e_results[0][1], "N/A (start)"))
    for i in range(1, 4):
        t1, p1 = high_e_results[i-1]
        t2, p2 = high_e_results[i]
        dt = t2 - t1
        rate = (p2 - p1) / dt
        print("{:6.1f} | {:14.2f} | {:<35.2f}".format(t2, p2, rate))
    print("The rate is decreasing rapidly, so there is no clear linear phase.\n")
    
    # --- 2. Simulate the corrected low enzyme concentration ---
    print("--- SCENARIO 2: Low Enzyme Concentration (E = {:.1f} uM) ---".format(E_low))
    print("Result: The reaction is slower, and the initial rate is constant (linear).")
    low_e_results = run_simulation(E_low, S_initial, Km, kcat, total_time, time_steps)
    
    # Print a few time points to show the linear curve
    print("Time (s) | Product (uM) | Rate (uM/s) in previous interval")
    print("-------- | -------------- | -----------------------------------")
    print("{:6.1f} | {:14.2f} | {}".format(low_e_results[0][0], low_e_results[0][1], "N/A (start)"))
    for i in range(1, 4):
        t1, p1 = low_e_results[i-1]
        t2, p2 = low_e_results[i]
        dt = t2 - t1
        rate = (p2 - p1) / dt
        print("{:6.1f} | {:14.2f} | {:<35.2f}".format(t2, p2, rate))
    print("The rate is nearly constant, showing a clear linear phase.\n")
    
    # --- 3. Demonstrate the calculation of the linear initial rate ---
    print("--- CALCULATING INITIAL VELOCITY ---")
    print("With a linear phase, we can now accurately calculate the initial velocity (V₀).")
    t_1, p_1 = low_e_results[1] # e.g., time point 1
    t_2, p_2 = low_e_results[2] # e.g., time point 2
    
    initial_velocity = (p_2 - p_1) / (t_2 - t_1)
    
    print("Using two points from the early, linear part of the low-enzyme reaction:")
    print(f"Point 1: (Time = {t_1:.1f} s, Product = {p_1:.2f} uM)")
    print(f"Point 2: (Time = {t_2:.1f} s, Product = {p_2:.2f} uM)")
    print(f"Initial Velocity = (Product_2 - Product_1) / (Time_2 - Time_1)")
    print(f"Initial Velocity = ({p_2:.2f} - {p_1:.2f}) / ({t_2:.1f} - {t_1:.1f})")
    print(f"Initial Velocity = {initial_velocity:.2f} uM/s")

if __name__ == "__main__":
    main()
