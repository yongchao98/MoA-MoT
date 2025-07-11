import numpy as np
from scipy.integrate import solve_ivp

def estimate_blowup_set_area():
    """
    This function numerically estimates the area of the set Omega by simulating
    the ODE system for a grid of initial conditions in the domain [-1, 1] x [2, 3].
    """
    
    # 1. Define the ODE system for a(t) and b(t)
    def ode_system(t, y):
        a, b = y
        da_dt = -b * a
        db_dt = -b**2 / 2 - a**2 + 6 * (a - 1)
        return [da_dt, db_dt]

    # 2. Set up the grid of initial conditions (a0, b0)
    # An odd number for `na` ensures that a0=0 is sampled.
    na = 101  # Number of points for a0 in [-1, 1]
    nb = 51   # Number of points for b0 in [2, 3]
    a0_range = np.linspace(-1, 1, na)
    b0_range = np.linspace(2, 3, nb)

    # 3. Initialize counters and simulation parameters
    t_span = [0, 5.0]  # Integration time window
    # We use a threshold to identify numerical blow-up
    blowup_threshold = 1e4
    
    blowup_count = 0
    total_points = na * nb

    # 4. Iterate over the grid, solve the ODE, and check for blow-up
    for a0 in a0_range:
        for b0 in b0_range:
            y0 = [a0, b0]
            
            # As reasoned in the analysis, we only expect blow-up to positive infinity if a(0) > 0.
            # We simulate all points to be thorough, but we only expect blow-up_count
            # to increment when a0 > 0.
            
            # Solve the ODE for the initial condition (a0, b0)
            sol = solve_ivp(
                ode_system,
                t_span,
                y0,
                method='RK45'
            )
            
            # Extract final values of a and b
            final_a = sol.y[0, -1]
            final_b = sol.y[1, -1]
            
            # Check if the trajectory meets the blow-up condition.
            # This happens if 'a' becomes large and positive, and 'b' becomes large and negative.
            # The solver might also fail (status < 0) near a singularity. We check for this too.
            if (sol.status < 0 and final_a > 0 and final_b < 0) or \
               (final_a > blowup_threshold and final_b < -blowup_threshold):
                blowup_count += 1
    
    # 5. Calculate the estimated measure of Omega
    # The total area of the initial sampling domain [-1, 1] x [2, 3]
    total_area = (a0_range[-1] - a0_range[0]) * (b0_range[-1] - b0_range[0])
    
    # The proportion of grid points that belong to Omega
    proportion = blowup_count / total_points
    
    # The estimated measure (area) of the set Omega
    estimated_measure = proportion * total_area
    
    # 6. Print the results of the calculation step-by-step
    print(f"Number of initial grid points causing blow-up (a -> +inf, b -> -inf): {blowup_count}")
    print(f"Total number of grid points simulated: {total_points}")
    print(f"Total area of the domain D = [-1,1] x [2,3] is: {total_area}")
    print("\nFinal calculation:")
    print(f"m(Omega) ~= (blow-up count / total points) * Total Area")
    print(f"m(Omega) ~= ({blowup_count} / {total_points}) * {total_area:.1f}")
    print(f"m(Omega) ~= {proportion:.4f} * {total_area:.1f}")
    print(f"Estimated m(Omega) = {estimated_measure:.4f}")

# Run the estimation function
estimate_blowup_set_area()