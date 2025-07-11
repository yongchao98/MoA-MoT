import numpy as np

def solve_point_kinetics():
    """
    Solves the point kinetics equations for a step reactivity insertion
    to demonstrate time evolution of reactor power.
    
    The equations are:
    dP/dt = (rho - beta) / Lambda * P + lambda_c * C
    dC/dt = beta / Lambda * P - lambda_c * C
    
    where:
    P = Power
    C = Delayed neutron precursor concentration
    rho = Reactivity
    beta = Total delayed neutron fraction
    Lambda = Prompt neutron generation time
    lambda_c = Effective precursor decay constant
    """
    
    # --- Reactor Physics Parameters ---
    # These are the "numbers in the final equation"
    rho = 0.001  # Reactivity insertion (e.g., control rod withdrawal), in absolute units
    beta = 0.0065  # Total delayed neutron fraction for U-235
    Lambda = 1e-5  # Prompt neutron generation time in seconds
    lambda_c = 0.08  # Effective one-group precursor decay constant in 1/seconds
    
    # --- Simulation Parameters ---
    P0 = 100.0  # Initial power in percent
    C0 = beta / (Lambda * lambda_c) * P0  # Initial precursor concentration at equilibrium
    dt = 0.01  # Time step in seconds
    total_time = 5.0  # Total simulation time in seconds
    n_steps = int(total_time / dt)

    # --- Print Parameters ---
    print("--- Point Kinetics Simulation Parameters ---")
    print(f"Reactivity (rho): {rho}")
    print(f"Delayed Neutron Fraction (beta): {beta}")
    print(f"Prompt Neutron Generation Time (Lambda): {Lambda} s")
    print(f"Precursor Decay Constant (lambda_c): {lambda_c} 1/s")
    print("-------------------------------------------\n")

    # --- Initial Conditions ---
    P = P0
    C = C0
    time = 0.0
    
    print("--- Time Evolution of Reactor Power ---")
    print(f"Time: {time:.2f} s, Power: {P:.2f} %")
    
    # --- Time Evolution using Euler method ---
    for i in range(1, n_steps + 1):
        # Calculate derivatives
        dP_dt = ((rho - beta) / Lambda) * P + lambda_c * C
        dC_dt = (beta / Lambda) * P - lambda_c * C
        
        # Update power and precursor concentration
        P += dP_dt * dt
        C += dC_dt * dt
        time += dt
        
        # Print output at selected intervals
        if (i * dt) % 1.0 == 0:
            print(f"Time: {time:.2f} s, Power: {P:.2f} %")

solve_point_kinetics()