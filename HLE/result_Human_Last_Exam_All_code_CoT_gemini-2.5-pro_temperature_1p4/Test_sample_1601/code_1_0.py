import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_and_estimate_area():
    """
    This function estimates the size of the set Omega based on numerical simulations.
    """
    # 1. Define the ODE system
    # y[0] = b(t), y[1] = a(t)
    def ode_system(t, y):
        b, a = y
        # Cap exp(t) to prevent overflow for large t, its qualitative effect remains.
        et = math.exp(t) if t < 100 else math.inf
        dbdt = -0.5 * b**2 - et * a**2 - a
        dadt = -b * a
        return [dbdt, dadt]

    # 2. Define the domain and total area
    a0_range = [-10.0, 1.0]
    b0_range = [10.0, 20.0]
    total_area = (a0_range[1] - a0_range[0]) * (b0_range[1] - b0_range[0])

    # 3. Create a grid of cell-center sample points for accurate area estimation
    # We choose n_a to be a multiple of the width of the a-range (11)
    # and n_b to be a multiple of the width of the b-range (10).
    n_a = 55
    n_b = 50

    a_step = (a0_range[1] - a0_range[0]) / n_a
    b_step = (b0_range[1] - b0_range[0]) / n_b

    a0_vals = np.linspace(a0_range[0] + a_step/2, a0_range[1] - a_step/2, n_a)
    b0_vals = np.linspace(b0_range[0] + b_step/2, b0_range[1] - b_step/2, n_b)

    # Define an event function to detect blow-up and terminate integration
    UPPER_BOUND = 1e6
    def event_stop(t, y):
        return max(abs(y[0]), abs(y[1])) - UPPER_BOUND
    event_stop.terminal = True

    t_span = [0, 10]
    blowup_count = 0
    total_count = n_a * n_b

    # 4. Loop through the grid of initial conditions
    for a0 in a0_vals:
        # 5. Skip non-candidate points based on analysis (a(t) can't go to +inf if a0 <= 0)
        if a0 <= 0:
            continue
        
        # Test all b0 values for this a0
        for b0 in b0_vals:
            y0 = [b0, a0]
            
            # 6. Solve the ODE for each candidate point
            sol = solve_ivp(
                ode_system,
                t_span,
                y0,
                events=event_stop,
                max_step=0.1
            )
            
            # 7. Check for confirmed blow-up condition
            # We check if the integration was terminated by our event.
            if sol.status == 1:
                final_b, final_a = sol.y[:, -1]
                # Confirm the trajectory is correct: a is large positive, b is negative
                if final_a > 0 and final_b < 0:
                    blowup_count += 1
    
    # 8. Estimate the area of Omega
    estimated_area = (blowup_count / total_count) * total_area

    # 9. Print the final calculation as requested
    print(f"Number of initial points leading to blow-up: {blowup_count}")
    print(f"Total number of initial points sampled: {total_count}")
    print(f"Total area of initial region: {total_area:.2f}")
    print(f"Estimated size of Omega, m(Omega):")
    print(f"m(Omega) = ({blowup_count} / {total_count}) * {total_area:.2f} = {estimated_area:.4f}")

solve_and_estimate_area()