import numpy as np
from scipy.integrate import solve_ivp

def solve_and_check_blowup():
    """
    This function implements the plan to numerically estimate the area of Omega.
    """
    
    # 1. Define the system of differential equations
    def system_of_odes(t, y):
        a, b = y
        a_prime = -b * a
        b_prime = -b**2 / 2 - a**2 + 6 * (a - 1)
        return [a_prime, b_prime]

    # 2. Set up simulation parameters and the grid of initial conditions
    a_range = [-1, 1]
    b_range = [2, 3]
    n_a = 101  # Use an odd number to include a=0 on the grid
    n_b = 51
    t_span = [0, 10]
    blowup_threshold = 1000.0  # Threshold to classify a solution as "blown up"

    a_vals = np.linspace(a_range[0], a_range[1], n_a)
    b_vals = np.linspace(b_range[0], b_range[1], n_b)
    
    total_points = n_a * n_b
    omega_count = 0

    # 3. Iterate through each initial condition and solve the ODE
    for a0 in a_vals:
        for b0 in b_vals:
            y0 = [a0, b0]
            
            # Solve the ODE. Use dense_output=False for efficiency as we only need the final point.
            # max_step helps handle the stiffness near blow-up without failing.
            sol = solve_ivp(system_of_odes, t_span, y0, max_step=0.1)
            
            # 4. Check the final state for blow-up condition
            a_final = sol.y[0, -1]
            b_final = sol.y[1, -1]
            
            if a_final > blowup_threshold and b_final < -blowup_threshold:
                omega_count += 1
    
    # 5. Estimate the area of Omega and print the results
    total_area = (a_range[1] - a_range[0]) * (b_range[1] - b_range[0])
    fraction_omega = omega_count / total_points
    estimated_area = fraction_omega * total_area
    
    print("Numerical Estimation Results:")
    print(f"Number of initial points in grid: {total_points}")
    print(f"Number of points in Omega (blow-up): {omega_count}")
    print(f"Fraction of Omega points: {fraction_omega:.4f}")
    print("---")
    print("Final Area Calculation:")
    print(f"Total area of initial region [-1,1]x[2,3]: {total_area}")
    print(f"Estimated m(Omega) = (Fraction of Omega points) * (Total Area)")
    print(f"m(Omega) = {fraction_omega:.4f} * {total_area}")
    print(f"Final estimated value for m(Omega) = {estimated_area:.4f}")

# Execute the simulation
solve_and_check_blowup()