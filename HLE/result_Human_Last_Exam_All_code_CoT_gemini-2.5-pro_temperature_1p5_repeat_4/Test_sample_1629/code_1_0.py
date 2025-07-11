import numpy as np
from scipy.integrate import solve_ivp

def estimate_blowup_set_area():
    """
    Estimates the area of the set of initial conditions leading to a specific
    blow-up behavior for a given system of ODEs.
    """

    # Define the system of ODEs
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # Set up the domain for initial conditions
    a_min, a_max = -1.0, 1.0
    b_min, b_max = 2.0, 3.0
    total_area = (a_max - a_min) * (b_max - b_min)

    # Set up the grid of initial points
    # Using an even number for num_a_points to avoid a point exactly at a=0,
    # simplifying the ratio calculation.
    num_a_points = 100
    num_b_points = 50
    a_points = np.linspace(a_min, a_max, num_a_points)
    b_points = np.linspace(b_min, b_max, num_b_points)
    total_points = num_a_points * num_b_points

    # Integration parameters
    t_span = [0, 20]  # Max integration time
    blowup_threshold = 100.0  # Threshold to detect blow-up

    # Event function to stop integration when solution gets large
    def event_blowup(t, y):
        a, b = y
        return max(abs(a), abs(b)) - blowup_threshold
    event_blowup.terminal = True  # Stop integration when event occurs

    # Loop through grid and count blow-up points
    blowup_counter = 0
    for a0 in a_points:
        for b0 in b_points:
            y0 = [a0, b0]
            # Solve the ODE
            sol = solve_ivp(
                ode_system,
                t_span,
                y0,
                events=event_blowup,
                method='RK45',
                dense_output=True
            )

            # Check if integration was terminated by the blow-up event
            if sol.status == 1 and sol.t[-1] < t_span[1]:
                # Get the final state
                a_final, b_final = sol.y[:, -1]
                # Check if it matches the specified blow-up condition
                if a_final > 0 and b_final < 0:
                    blowup_counter += 1

    # Estimate the area of the blow-up set Omega
    estimated_area = (blowup_counter / total_points) * total_area
    
    print(f"Numerical Estimation of m(Omega)")
    print("---------------------------------")
    print(f"Domain of initial conditions: a in [{a_min}, {a_max}], b in [{b_min}, {b_max}]")
    print(f"Total area of domain: {total_area}")
    print(f"Grid size: {num_a_points} x {num_b_points} = {total_points} total points")
    print(f"Number of initial points leading to blow-up: {blowup_counter}")
    print("\nFinal Area Calculation:")
    print(f"m(Omega) ~= (blowup_points / total_points) * total_area")
    print(f"m(Omega) ~= ({blowup_counter} / {total_points}) * {total_area}")
    print(f"Estimated m(Omega) = {estimated_area}")


if __name__ == '__main__':
    estimate_blowup_set_area()
