import numpy as np

def simulate_enzyme_kinetics(E_conc, S_initial, K_m, k_cat, total_time, time_step):
    """
    Simulates product formation over time in an enzyme kinetics assay.

    Args:
        E_conc (float): Enzyme concentration in µM.
        S_initial (float): Initial substrate concentration in µM.
        K_m (float): Michaelis constant in µM.
        k_cat (float): Catalytic rate constant (turnover number) in 1/s.
        total_time (int): Total simulation time in seconds.
        time_step (float): Time step for the simulation in seconds.

    Returns:
        tuple: A tuple containing lists of time points, product concentrations,
               and substrate concentrations.
    """
    # Calculate Vmax from enzyme concentration and kcat
    V_max = k_cat * E_conc

    print(f"\n--- Simulation with [E] = {E_conc:.2f} µM ---")
    # Using the equation V_max = k_cat * [E]
    print(f"Equation parameters: V_max = {k_cat} s⁻¹ * {E_conc:.2f} µM = {V_max:.2f} µM/s, K_m = {K_m} µM")

    # Initialize lists to store data
    time_points = [0]
    substrate_conc = [S_initial]
    product_conc = [0]

    current_S = S_initial
    current_P = 0

    for t in np.arange(time_step, total_time + time_step, time_step):
        # Calculate reaction velocity using the Michaelis-Menten equation
        # v = (V_max * [S]) / (K_m + [S])
        velocity = (V_max * current_S) / (K_m + current_S)

        # Calculate substrate consumed and product formed in this time step
        product_formed_in_step = velocity * time_step
        
        # Update concentrations
        current_S -= product_formed_in_step
        # Ensure substrate doesn't go below zero
        current_S = max(0, current_S)
        current_P += product_formed_in_step

        time_points.append(round(t, 2))
        substrate_conc.append(round(current_S, 2))
        product_conc.append(round(current_P, 2))

    return time_points, product_conc

def main():
    """Main function to run simulations and print results."""
    # --- Assay Parameters ---
    S_INITIAL = 100.0  # µM
    K_M = 20.0       # µM
    K_CAT = 10.0       # 1/s
    TOTAL_TIME = 60    # seconds
    TIME_STEP = 5      # seconds

    print("Problem: The 'Product vs Time' plot is non-linear, suggesting the reaction is too fast.")
    print("This simulation will demonstrate why decreasing enzyme concentration helps.")

    # --- Scenario 1: High Enzyme Concentration (The Problem) ---
    high_E_conc = 1.0  # µM
    time1, product1 = simulate_enzyme_kinetics(high_E_conc, S_INITIAL, K_M, K_CAT, TOTAL_TIME, TIME_STEP)
    
    print("Product vs. Time data ([Product] in µM):")
    for t, p in zip(time1, product1):
        print(f"Time: {t:2d} s, Product: {p:5.2f} µM")
    print("Observation: The rate of product formation decreases significantly over time. The plot is not linear.")

    # --- Scenario 2: Low Enzyme Concentration (The Solution) ---
    low_E_conc = 0.1   # µM
    time2, product2 = simulate_enzyme_kinetics(low_E_conc, S_INITIAL, K_M, K_CAT, TOTAL_TIME, TIME_STEP)
    
    print("Product vs. Time data ([Product] in µM):")
    for t, p in zip(time2, product2):
        print(f"Time: {t:2d} s, Product: {p:5.2f} µM")
    print("Observation: The increase in product is much more constant over time.")
    print("The plot is now linear, which is ideal for measuring the initial reaction rate.")

if __name__ == "__main__":
    main()