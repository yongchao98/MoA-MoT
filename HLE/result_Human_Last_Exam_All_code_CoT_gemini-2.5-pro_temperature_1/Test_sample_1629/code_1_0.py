import numpy as np
from scipy.integrate import solve_ivp

def estimate_blowup_set_measure():
    """
    Numerically estimates the measure of the set Omega.
    """
    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Define an event to detect blow-up, stopping the integration early.
    # The event is triggered when |a| or |b| exceeds a large threshold.
    blowup_threshold = 1000.0
    def blowup_event(t, y):
        return max(abs(y[0]), abs(y[1])) - blowup_threshold
    blowup_event.terminal = True  # Stop integration when the event occurs

    # Set up the grid of initial conditions in the domain [-1, 1] x [2, 3]
    grid_resolution = 200  # Use a 200x200 grid for good accuracy
    a0_vals = np.linspace(-1, 1, grid_resolution)
    b0_vals = np.linspace(2, 3, grid_resolution)

    blowup_count = 0
    total_count = grid_resolution * grid_resolution

    # Iterate over each initial condition in the grid
    for a0 in a0_vals:
        for b0 in b0_vals:
            initial_condition = [a0, b0]
            
            # Solve the ODE for this initial condition
            solution = solve_ivp(
                ode_system, 
                [0, 20],  # Max integration time
                initial_condition, 
                method='RK45', 
                events=blowup_event,
                dense_output=True
            )
            
            # Check if the integration was terminated by the blow-up event
            if solution.status == 1:
                # Get the state at the moment the event was triggered
                a_final, b_final = solution.y_events[0][-1]
                
                # Check if the blow-up matches the condition for Omega
                # (a -> +inf, b -> -inf)
                if a_final > (blowup_threshold - 1) and b_final < -(blowup_threshold - 1):
                    blowup_count += 1

    # The area of the initial domain R = [-1, 1] x [2, 3]
    domain_width = 1 - (-1)
    domain_height = 3 - 2
    domain_area = domain_width * domain_height
    
    # Estimate the measure of Omega
    fraction_in_omega = blowup_count / total_count
    m_omega_estimate = fraction_in_omega * domain_area
    
    print("Numerical Estimation of m(Omega):")
    print(f"Domain R = [-1, 1] x [2, 3], Area = {domain_width} * {domain_height} = {domain_area}")
    print(f"Grid size: {grid_resolution}x{grid_resolution} = {total_count} points")
    print(f"Initial points leading to specified blow-up: {blowup_count}")
    print(f"Fraction of points in Omega: {blowup_count} / {total_count} = {fraction_in_omega}")
    print(f"Estimated measure m(Omega) = fraction * Area = {fraction_in_omega} * {domain_area} = {m_omega_estimate}")

estimate_blowup_set_measure()