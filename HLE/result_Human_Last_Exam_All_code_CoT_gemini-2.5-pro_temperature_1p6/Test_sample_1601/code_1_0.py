import numpy as np
from scipy.integrate import solve_ivp

def solve_system_and_check_blowup():
    """
    This function solves the given system of ODEs for a grid of initial conditions
    and estimates the measure of the set Omega that leads to a specific blow-up behavior.
    """

    # Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - np.exp(t) * a**2 - a
        return [da_dt, db_dt]

    # Define an event to stop integration if a variable blows up
    # This speeds up the simulation by not integrating past the blow-up.
    BLOWUP_THRESHOLD = 1e6
    def blowup_event(t, y):
        return np.max(np.abs(y)) - BLOWUP_THRESHOLD
    blowup_event.terminal = True
    blowup_event.direction = 1

    # Define the domain for initial conditions
    a0_min, a0_max = -10.0, 1.0
    b0_min, b0_max = 10.0, 20.0
    domain_area = (a0_max - a0_min) * (b0_max - b0_min)

    # Set up a representative grid of initial conditions using cell-centers
    # This ensures the proportion of test points in a region matches the proportion of its area.
    na = 110  # Number of points for a0
    nb = 100  # Number of points for b0
    
    da = (a0_max - a0_min) / na
    db = (b0_max - b0_min) / nb
    
    a0_vals = np.linspace(a0_min, a0_max, na, endpoint=False) + da / 2
    b0_vals = np.linspace(b0_min, b0_max, nb, endpoint=False) + db / 2
    
    total_points = na * nb
    successful_count = 0
    
    t_span = [0, 10] # Time interval for integration

    # Loop through all initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            # Initial condition vector
            y0 = [a0, b0]

            # Solve the ODE
            sol = solve_ivp(ode_system, t_span, y0, events=blowup_event, dense_output=True)
            
            # Check the status of the solution at the end of integration
            a_final, b_final = sol.y_at(sol.t_events[0])[0] if sol.t_events[0].size > 0 else sol.y[:,-1]

            # Check if the solution matches the specified blow-up condition
            if a_final > BLOWUP_THRESHOLD / 2 and b_final < -BLOWUP_THRESHOLD / 2:
                successful_count += 1
    
    # Estimate the measure of Omega
    if total_points > 0:
        estimated_area = (successful_count / total_points) * domain_area
    else:
        estimated_area = 0

    # Print the result in the required format
    print(f"Numerical Estimation of m(Omega):")
    print(f"Total initial conditions sampled: {total_points}")
    print(f"Initial conditions leading to blow-up: {successful_count}")
    print(f"Total area of the initial domain: {domain_area}")
    print(f"The equation for the estimated measure is:")
    print(f"m(Omega) = ({successful_count} / {total_points}) * {domain_area}")
    print(f"Result: m(Omega) ~= {estimated_area}")

if __name__ == '__main__':
    solve_system_and_check_blowup()