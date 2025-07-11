import numpy as np
from scipy.integrate import solve_ivp

def estimate_blowup_set_measure():
    """
    Estimates the measure of the set Omega by numerically integrating the ODE
    system for a grid of initial conditions and counting the fraction that
    leads to the specified blow-up.
    """

    # 1. Define the ODE system
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # 2. Define the problem domain and parameters
    a_lims = [-1, 1]
    b_lims = [2, 3]
    domain_area = (a_lims[1] - a_lims[0]) * (b_lims[1] - b_lims[0])

    # 3. Set up the grid of initial conditions
    N_a = 100  # Number of points along the a-axis
    N_b = 50   # Number of points along the b-axis
    total_points = N_a * N_b
    a0_vals = np.linspace(a_lims[0], a_lims[1], N_a)
    b0_vals = np.linspace(b_lims[0], b_lims[1], N_b)

    # 4. Set up integration and event detection parameters
    t_max = 20.0
    blowup_threshold = 100.0

    def escape_event(t, y):
        # Event function is zero when the solution hits the boundary
        return max(abs(y[0]), abs(y[1])) - blowup_threshold
    escape_event.terminal = True  # Stop integration when event occurs

    # 5. Loop over initial conditions and solve
    blowup_count = 0
    for a0 in a0_vals:
        # The line a=0 is a set of measure zero, we can skip it.
        # Trajectories starting with a(0)<0 are expected to result in a -> -inf.
        if a0 <= 0:
            continue
            
        for b0 in b0_vals:
            y0 = [a0, b0]
            
            sol = solve_ivp(
                ode_system,
                [0, t_max],
                y0,
                events=escape_event,
                rtol=1e-6, atol=1e-6
            )
            
            # Check if the integration was stopped by the event
            if sol.status == 1:
                # Get the state vector at the time of the event
                a_final, b_final = sol.y_events[0][:, 0]
                
                # Check if it escaped into the quadrant where a>0 and b<0
                if a_final > 0 and b_final < 0:
                    blowup_count += 1
    
    # 6. Estimate the area m(Omega)
    # The domain for which we expect blowup is a>0, which is half the grid.
    expected_blowup_points = (N_a // 2) * N_b
    
    # The area is estimated by the fraction of points that led to blow-up
    # multiplied by the total area of the domain.
    # Note: Since we only tested points with a0 > 0, we adjust the calculation.
    # The area of the tested sub-domain (a>0) is 1.
    area_of_tested_subdomain = (1 - 0) * (b_lims[1] - b_lims[0])
    
    # The number of points in the tested sub-domain (a > 0).
    # np.linspace includes endpoints, so we find how many a0 are > 0.
    num_a_positive = np.sum(a0_vals > 0)
    total_positive_a_points = num_a_positive * N_b

    if total_positive_a_points > 0:
      m_omega_est = area_of_tested_subdomain * (blowup_count / total_positive_a_points)
    else:
      m_omega_est = 0


    # 7. Print the results clearly
    print("Numerical Estimation of m(Omega)")
    print("---------------------------------")
    print(f"The simulation tested initial conditions in the region (0, 1] x [2, 3].")
    print(f"Area of this sub-region is {area_of_tested_subdomain}.")
    print(f"Number of grid points in this sub-region: {total_positive_a_points}")
    print(f"Number of initial conditions leading to blow-up (a -> +inf, b -> -inf): {blowup_count}")
    print("\nFinal Calculation:")
    print(f"m(Omega) ~= (Area of Sub-Region) * (blow-up points / total points in sub-region)")
    print(f"m(Omega) ~= {area_of_tested_subdomain} * ({blowup_count} / {total_positive_a_points}) = {m_omega_est:.4f}")
    print(f"\nThe estimated measure m(Omega) is approximately {round(m_omega_est)}.")

if __name__ == '__main__':
    estimate_blowup_set_measure()