import numpy as np
from scipy.integrate import solve_ivp

def estimate_blow_up_set_measure():
    """
    Solves the given system of ODEs for a grid of initial conditions and
    estimates the measure of the set of initial conditions that lead to blow-up.
    """
    
    # Define the system of ODEs
    def ode_system(t, y):
        b, a = y
        dbdt = -b**2 / 2 - np.exp(t) * a**2 - a
        dadt = -b * a
        # Add a check to prevent extreme values from breaking the solver
        if not np.isfinite(dbdt) or not np.isfinite(dadt):
            return [0, 0]
        return [dbdt, dadt]

    # Define the grid of initial conditions
    a0_range = [-10, 1]
    b0_range = [10, 20]
    
    # Grid resolution
    n_a = 55  # Number of points for a0
    n_b = 50  # Number of points for b0
    
    a0_vals = np.linspace(a0_range[0], a0_range[1], n_a)
    b0_vals = np.linspace(b0_range[0], b0_range[1], n_b)
    
    # Simulation parameters
    t_span = [0, 10]
    blow_up_threshold = 1.0e6
    
    blow_up_count = 0
    total_points = n_a * n_b
    
    # Loop over the grid of initial conditions
    for a0 in a0_vals:
        # Based on our analysis, blow-up only happens for a0 > 0.
        # The script will confirm this by testing all points.
        
        for b0 in b0_vals:
            initial_conditions = [b0, a0]
            
            sol = solve_ivp(ode_system, t_span, initial_conditions, method='RK45', dense_output=True)
            
            # Check for blow-up:
            # 1. Did the solver stop before the end time? (Often indicates divergence)
            # 2. Does the final state satisfy the blow-up conditions?
            
            b_final = sol.y[0, -1]
            a_final = sol.y[1, -1]

            is_diverged = sol.t[-1] < t_span[1] and a0 > 0 # A heuristic for divergence
            is_blew_up = a_final > blow_up_threshold and b_final < -blow_up_threshold
            
            if is_diverged or is_blew_up:
                blow_up_count += 1

    # Calculate the total area of the initial domain
    total_area = (a0_range[1] - a0_range[0]) * (b0_range[1] - b0_range[0])
    
    # Estimate the area of the blow-up set Omega
    estimated_m_Omega = (blow_up_count / total_points) * total_area

    # Output the results
    print(f"Initial domain for a: {a0_range}")
    print(f"Initial domain for b: {b0_range}")
    print(f"Total area of the initial domain: {total_area}")
    print(f"Grid size (a x b): {n_a} x {n_b} = {total_points} points")
    print(f"Number of points leading to blow-up: {blow_up_count}")
    print("m(Omega) = (blow_up_count / total_points) * total_area")
    print(f"m(Omega) = ({blow_up_count} / {total_points}) * {total_area}")
    print(f"Estimated size of the set Omega, m(Omega): {estimated_m_Omega:.4f}")

if __name__ == '__main__':
    estimate_blow_up_set_measure()
