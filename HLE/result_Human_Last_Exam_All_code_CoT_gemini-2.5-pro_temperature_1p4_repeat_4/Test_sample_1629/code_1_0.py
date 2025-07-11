import numpy as np
from scipy.integrate import solve_ivp

def estimate_omega_measure():
    """
    Estimates the measure of the set Omega by numerically integrating the ODE
    for a grid of initial conditions.
    """

    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Set up the grid of initial conditions in [-1, 1] x [2, 3]
    # Use an odd number for N_a to ensure a=0 is sampled
    N_a = 201 
    N_b = 101
    a0_vals = np.linspace(-1, 1, N_a)
    b0_vals = np.linspace(2, 3, N_b)

    # Integration parameters
    t_max = 50.0
    t_span = [0, t_max]

    # Blow-up condition thresholds
    A_LARGE = 1.0e4
    B_SMALL = -1.0e4

    # Counter for blow-up trajectories
    blow_up_count = 0
    total_points = N_a * N_b

    # Loop over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            y0 = [a0, b0]
            
            # Solve the ODE for the current initial condition
            sol = solve_ivp(
                ode_system, 
                t_span, 
                y0, 
                method='RK45', 
                dense_output=True,
                rtol=1e-6,
                atol=1e-6
            )
            
            # Check if the final state satisfies the blow-up condition
            final_a, final_b = sol.y[0, -1], sol.y[1, -1]
            if final_a > A_LARGE and final_b < B_SMALL:
                blow_up_count += 1
    
    # Calculate the area of the initial rectangle
    # area = (x_max - x_min) * (y_max - y_min)
    total_area_x_min = -1
    total_area_x_max = 1
    total_area_y_min = 2
    total_area_y_max = 3
    total_area = (total_area_x_max - total_area_x_min) * (total_area_y_max - total_area_y_min)

    # Estimate the measure of Omega
    # m(Omega) = (fraction of blow-ups) * (total area)
    area_of_omega = (blow_up_count / total_points) * total_area

    # Output the details of the final calculation
    print(f"Numerical Estimation Details:")
    print(f"Initial conditions grid: {N_a} x {N_b} = {total_points} points")
    print(f"Total area of initial rectangle: ({total_area_x_max} - ({total_area_x_min})) * ({total_area_y_max} - {total_area_y_min}) = {total_area}")
    print(f"Initial conditions leading to specified blow-up: {blow_up_count}")
    print(f"\nFinal Equation for Estimation:")
    print(f"m(Omega) = (blow_up_count / total_points) * total_area")
    print(f"m(Omega) = ({blow_up_count} / {total_points}) * {total_area}")
    print(f"\nEstimated size of Omega: {area_of_omega:.4f}")

if __name__ == '__main__':
    estimate_omega_measure()