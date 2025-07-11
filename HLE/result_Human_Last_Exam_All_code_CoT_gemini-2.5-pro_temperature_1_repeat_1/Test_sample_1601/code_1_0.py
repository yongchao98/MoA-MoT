import numpy as np
from scipy.integrate import solve_ivp

def solve_ode_and_estimate_area():
    """
    Solves the given system of ODEs for a grid of initial conditions
    and estimates the area of the set Omega where solutions blow up.
    """

    # Define the system of ordinary differential equations
    # y[0] corresponds to a(t), y[1] corresponds to b(t)
    def ode_system(t, y):
        a, b = y
        # The exponential term e^t can grow very large.
        # The solver will typically stop due to blow-up or event detection
        # before t gets large enough for this to be a numerical issue.
        da_dt = -b * a
        db_dt = -0.5 * b**2 - np.exp(t) * a**2 - a
        return [da_dt, db_dt]

    # Define an event function to detect when b(t) crosses zero.
    # The integration will terminate when this event occurs.
    def event_b_crosses_zero(t, y):
        return y[1]  # The event is when b(t) = 0
    event_b_crosses_zero.terminal = True  # Stop integration when event occurs
    event_b_crosses_zero.direction = -1  # Trigger only when b(t) goes from positive to negative

    # Set up the grid of initial conditions for a(0) and b(0).
    # We only need to check a(0) > 0.
    # A 50x50 grid provides a good balance of accuracy and computation time.
    Na = 50  # Number of points for a(0)
    Nb = 50  # Number of points for b(0)
    a0_values = np.linspace(0.02, 1.0, Na)  # a(0) in (0, 1]
    b0_values = np.linspace(10.0, 20.0, Nb) # b(0) in [10, 20]

    blowup_count = 0
    total_count = Na * Nb
    
    # Maximum integration time. If b(t) doesn't cross zero by this time,
    # it's assumed it never will, as a(t) will have decayed significantly.
    t_max = 10.0

    # Iterate over the grid of initial conditions
    for a0 in a0_values:
        for b0 in b0_values:
            initial_conditions = [a0, b0]
            
            # Numerically solve the ODE system
            solution = solve_ivp(
                ode_system,
                [0, t_max],
                initial_conditions,
                events=event_b_crosses_zero,
                dense_output=True # Required for accurate event location
            )
            
            # Check if the integration was terminated by our event
            # sol.status == 1 indicates a termination event occurred.
            if solution.status == 1 and solution.t_events[0].size > 0:
                blowup_count += 1
    
    # The area of the domain where blow-up is possible is (1-0) * (20-10)
    total_area = 10.0
    
    # Estimate the area of Omega
    omega_area_estimate = (blowup_count / total_count) * total_area

    print("--- Estimation of m(Omega) ---")
    print(f"Total initial conditions sampled: {total_count}")
    print(f"Initial conditions leading to blow-up: {blowup_count}")
    print(f"Total area of the search domain (0, 1] x [10, 20]: {total_area}")
    print("\nFinal area estimation equation:")
    # The problem asks to "output each number in the final equation"
    print(f"m(Omega) ~= ({blowup_count} / {total_count}) * {total_area}")
    
    print(f"\nEstimated size m(Omega): {omega_area_estimate:.4f}")

# Execute the function
solve_ode_and_estimate_area()