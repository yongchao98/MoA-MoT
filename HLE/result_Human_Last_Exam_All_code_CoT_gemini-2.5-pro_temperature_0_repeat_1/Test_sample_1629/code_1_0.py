import numpy as np
from scipy.integrate import solve_ivp

def solve_and_check():
    """
    This function performs the numerical estimation of the area of Omega.
    """
    # Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Event function for a(t) = 0. Terminal event.
    def event_a_zero(t, y):
        return y[0]
    event_a_zero.terminal = True

    # Event function for b(t) = 0. Terminal event.
    def event_b_zero(t, y):
        return y[1]
    event_b_zero.terminal = True

    # Set up the grid of initial conditions in the relevant domain (0, 1] x [2, 3]
    # The area of this domain is 1.0
    N = 200  # Grid resolution
    a0_vals = np.linspace(1e-6, 1.0, N)  # a0 > 0
    b0_vals = np.linspace(2.0, 3.0, N)

    omega_count = 0
    total_count = N * N

    # Time span for integration. Should be long enough for an event to occur.
    t_span = [0, 50]

    # Loop over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            y0 = [a0, b0]
            
            # Solve the ODE for the initial condition y0
            sol = solve_ivp(
                ode_system,
                t_span,
                y0,
                events=(event_a_zero, event_b_zero),
                method='RK45'
            )
            
            # Check if the integration was terminated by an event
            if sol.status == 1:
                t_a_events = sol.t_events[0]
                t_b_events = sol.t_events[1]
                
                # If b hits zero first, the point is in Omega.
                # This happens if t_b_events is not empty AND
                # (t_a_events is empty OR the first b event time is less than the first a event time)
                if len(t_b_events) > 0:
                    if len(t_a_events) == 0 or t_b_events[0] < t_a_events[0]:
                        omega_count += 1
    
    # The total area of the original domain is 2.
    # The area of the tested sub-domain (where Omega can exist) is 1.
    # m(Omega) is the fraction of the sub-domain that belongs to Omega.
    area_of_tested_domain = 1.0
    estimated_area = (omega_count / total_count) * area_of_tested_domain
    
    print("The analysis shows that solutions can only blow up as required if a(0) > 0.")
    print("We test the domain (0, 1] x [2, 3], which has an area of 1.0.")
    print(f"Number of initial points leading to blow-up: {omega_count}")
    print(f"Total number of initial points tested: {total_count}")
    print(f"Area of tested domain: {area_of_tested_domain}")
    print(f"Estimated measure m(Omega) = {omega_count} / {total_count} * {area_of_tested_domain} = {estimated_area}")

solve_and_check()