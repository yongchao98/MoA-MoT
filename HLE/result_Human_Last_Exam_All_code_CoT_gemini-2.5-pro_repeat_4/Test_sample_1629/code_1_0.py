import numpy as np
from scipy.integrate import solve_ivp

def estimate_omega_measure():
    """
    Estimates the measure of the set Omega for the given system of ODEs.
    """
    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -0.5 * b**2 - a**2 + 6 * (a - 1)
        return [dadt, dbdt]

    # Define the domain D for initial conditions (a0, b0)
    a0_range = [-1, 1]
    b0_range = [2, 3]
    domain_area = (a0_range[1] - a0_range[0]) * (b0_range[1] - b0_range[0])

    # Setup the grid for numerical simulation
    # A finer grid gives a more accurate estimate.
    num_a_points = 200
    num_b_points = 100
    total_points = num_a_points * num_b_points
    a0_vals = np.linspace(a0_range[0], a0_range[1], num_a_points)
    b0_vals = np.linspace(b0_range[0], b0_range[1], num_b_points)

    # Simulation parameters
    t_max = 10.0  # Time horizon for integration
    blowup_threshold = 1e6  # Threshold to detect numerical blow-up

    # Counter for initial conditions in Omega
    omega_count = 0

    # Iterate over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            # Set initial condition
            y0 = [a0, b0]

            # Solve the ODE system
            sol = solve_ivp(ode_system, [0, t_max], y0, t_eval=[t_max])

            # Check if the solution blew up as specified
            a_final, b_final = sol.y[:, -1]
            if a_final > blowup_threshold and b_final < -blowup_threshold:
                omega_count += 1

    # Estimate the measure of Omega
    estimated_measure = (omega_count / total_points) * domain_area
    
    # Print the details of the calculation
    print(f"Domain of initial conditions: a in [{a0_range[0]}, {a0_range[1]}], b in [{b0_range[0]}, {b0_range[1]}]")
    print(f"Total area of the domain: {domain_area}")
    print(f"Grid size: {num_a_points} x {num_b_points} = {total_points} total points")
    print(f"Number of points in Omega (leading to blow-up): {omega_count}")
    print("\nFinal calculation:")
    print(f"m(Omega) ~= (Number of points in Omega / Total points) * Total Area")
    print(f"m(Omega) ~= ({omega_count} / {total_points}) * {domain_area}")
    print(f"m(Omega) ~= {estimated_measure}")

if __name__ == '__main__':
    estimate_omega_measure()