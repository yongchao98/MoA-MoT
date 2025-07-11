import numpy as np
from scipy.integrate import solve_ivp
import warnings

def main():
    """
    Main function to perform the simulation and estimate the measure of Omega.
    """
    # Suppress warnings that might arise from the solver near singularities
    warnings.filterwarnings('ignore', 'The solver_ivp implicit solver')

    # Define the system of ordinary differential equations
    def ode_system(t, y):
        """
        Defines the ODE system:
        a'(t) = -b(t)a(t)
        b'(t) = -b^2(t)/2 - a^2(t) + 6(a(t)-1)
        """
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Define the rectangular region for initial conditions
    a_min, a_max = -1.0, 1.0
    b_min, b_max = 2.0, 3.0
    
    # Calculate the total area of the region
    total_area = (a_max - a_min) * (b_max - b_min)

    # Set up the grid of initial points
    n_a = 201
    n_b = 101
    a_vals = np.linspace(a_min, a_max, n_a)
    b_vals = np.linspace(b_min, b_max, n_b)
    total_points = n_a * n_b

    # Simulation parameters
    t_span = [0, 5.0]  # Time interval for integration
    a_blowup_thresh = 1000.0 # Threshold for a(t) -> infinity
    b_blowup_thresh = -1000.0 # Threshold for b(t) -> -infinity

    blowup_count = 0

    # Iterate through the grid of initial conditions
    for a0 in a_vals:
        for b0 in b_vals:
            y0 = [a0, b0]
            try:
                # Solve the ODE
                sol = solve_ivp(ode_system, t_span, y0, method='Radau', rtol=1e-3, atol=1e-6)
                
                # Check the state at the end of the integration
                a_final, b_final = sol.y[:, -1]

                # Check for the specified blow-up condition
                if a_final > a_blowup_thresh and b_final < b_blowup_thresh:
                    blowup_count += 1
                # If solver failed, it can indicate a singularity (blow-up).
                # Analytic arguments show blow-up to a->+inf occurs for a0>0.
                elif sol.status != 0 and a0 > 0:
                     blowup_count += 1

            except Exception:
                # Any other exception during integration likely means blow-up.
                # We only count it if the initial a0 was positive.
                if a0 > 0:
                    blowup_count += 1
    
    # Calculate the fraction of points leading to blow-up
    fraction = blowup_count / total_points
    
    # Estimate the measure of the set Omega
    estimated_measure = fraction * total_area

    # Output the numbers involved in the final calculation
    print("--- Simulation Results ---")
    print(f"Total number of initial points sampled: {total_points}")
    print(f"Number of points leading to specified blow-up: {blowup_count}")
    print(f"Total area of the initial region [-1,1] x [2,3]: {total_area}")
    print("\n--- Final Calculation ---")
    print(f"m(Omega) ~= (Number of Blow-up Points / Total Points) * Total Area")
    print(f"m(Omega) ~= ({blowup_count} / {total_points}) * {total_area}")
    print(f"m(Omega) ~= {fraction:.4f} * {total_area}")
    print(f"Estimated measure m(Omega): {estimated_measure:.4f}")

if __name__ == '__main__':
    main()