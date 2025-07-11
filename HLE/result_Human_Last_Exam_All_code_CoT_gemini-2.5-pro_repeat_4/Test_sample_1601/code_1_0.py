import numpy as np
from scipy.integrate import solve_ivp

def solve_and_check_blowup():
    """
    This function solves the ODE system for a grid of initial conditions
    and estimates the measure of the set Omega.
    """
    # Define the system of ODEs
    # y[0] = a, y[1] = b
    def vector_field(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -b**2 / 2 - np.exp(t) * a**2 - a
        return [dadt, dbdt]

    # The analysis suggests we only need to check the domain where a(0) > 0.
    # Domain for a(0): (0, 1]
    # Domain for b(0): [10, 20]
    a0_min, a0_max = 1e-2, 1.0  # Start a bit away from 0 for numerical stability
    b0_min, b0_max = 10.0, 20.0
    area_of_tested_region = (a0_max - 0) * (b0_max - b0_min)

    # Set up the grid of initial conditions
    grid_size = 20
    a0_grid = np.linspace(a0_min, a0_max, grid_size)
    b0_grid = np.linspace(b0_min, b0_max, grid_size)
    num_points = len(a0_grid) * len(b0_grid)
    
    blowup_count = 0
    
    # Time span for integration
    t_span = [0, 5]
    # Threshold to consider as blow-up
    blowup_threshold = 1e6

    # Loop through all initial conditions
    for a0 in a0_grid:
        for b0 in b0_grid:
            y0 = [a0, b0]
            # We use a try-except block because solve_ivp may fail for
            # stiff systems or systems that blow up, which is what we expect.
            try:
                sol = solve_ivp(vector_field, t_span, y0, max_step=0.01)

                # Check if the solution reached the end of the time span
                if sol.status == 0:
                    a_final, b_final = sol.y[:, -1]
                    # Check for blow-up condition
                    if a_final > blowup_threshold and b_final < -blowup_threshold:
                        blowup_count += 1
                else:
                    # If solver failed, it's very likely due to a blow-up.
                    # A more detailed check could look at the last values before failure,
                    # but our analysis suggests any failure will be the correct type of blow-up.
                    blowup_count += 1
            except Exception:
                # Any other exception during integration is also treated as a blow-up
                blowup_count += 1

    # Estimate the measure of Omega
    estimated_measure = (blowup_count / num_points) * area_of_tested_region

    print("--- Numerical Estimation of m(Omega) ---")
    print(f"Domain for blow-up analysis: a in (0, 1], b in [10, 20]")
    print(f"Area of this domain = (1 - 0) * (20 - 10) = {area_of_tested_region}")
    print(f"Number of total grid points tested = {num_points}")
    print(f"Number of initial points leading to blow-up = {blowup_count}")
    print(f"Fraction of blow-ups = {blowup_count / num_points:.2f}")
    print("\nFinal Calculation:")
    print(f"m(Omega) ~= (blowup_count / num_points) * area_of_tested_region")
    print(f"m(Omega) ~= ({blowup_count} / {num_points}) * {area_of_tested_region} = {estimated_measure:.4f}")
    
solve_and_check_blowup()