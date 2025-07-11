import numpy as np
from scipy.integrate import solve_ivp

def main():
    """
    This script estimates the measure of the set Omega for the given system of ODEs.
    """
    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -0.5 * b**2 - a**2 + 6 * a - 6
        return [dadt, dbdt]

    # 1. Define the domain of initial conditions
    a_range = [-1.0, 1.0]
    b_range = [2.0, 3.0]
    domain_area = (a_range[1] - a_range[0]) * (b_range[1] - b_range[0])

    # 2. Create a grid of initial conditions
    # Using a reasonably fine grid for accuracy. n_a is odd to include a0=0.
    n_a = 101
    n_b = 51
    a0_values = np.linspace(a_range[0], a_range[1], n_a)
    b0_values = np.linspace(b_range[0], b_range[1], n_b)
    total_points = n_a * n_b

    # 3. & 4. Loop through the grid, solve ODE, and check for blow-up condition
    t_max = 5.0  # A sufficiently long time for blow-up to occur
    omega_count = 0

    for a0 in a0_values:
        for b0 in b0_values:
            initial_conditions = [a0, b0]
            
            # Solve the ODE. 'RK45' is a good general-purpose solver that can
            # detect stiffness and singularities, often stopping before the end of the time interval.
            solution = solve_ivp(
                ode_system,
                [0, t_max],
                initial_conditions,
                method='RK45'
            )

            # Check for finite-time blow-up.
            # This is indicated if the solver did not reach the end of the time interval.
            # sol.status == 1 indicates that a singularity was encountered.
            if solution.status == 1:
                # Get the state at the last computed time step
                a_final = solution.y[0, -1]
                b_final = solution.y[1, -1]

                # The desired blow-up is a -> +inf and b -> -inf.
                # This means at the point of blow-up, 'a' must be positive and 'b' must be negative.
                if a_final > 0 and b_final < 0:
                    omega_count += 1
    
    # 5. Calculate and print the final result
    estimated_measure = (omega_count / total_points) * domain_area

    print("--- Estimation of the measure of set Omega ---")
    print(f"Domain for initial conditions (a(0), b(0)): {a_range} x {b_range}")
    print(f"Total area of the domain: {domain_area}")
    print(f"Grid size for numerical integration: {n_a}x{n_b} = {total_points} points")
    print("\nCondition for set Omega: a(t) -> +inf and b(t) -> -inf in finite time.")
    print(f"Number of initial points found in Omega: {omega_count}")
    print(f"Total number of initial points sampled: {total_points}")
    print("\nFinal calculation for the measure of Omega, m(Omega):")
    print(f"m(Omega) ≈ (Points in Omega / Total Points) * Domain Area")
    # Output each number in the final equation as requested
    print(f"m(Omega) ≈ ({omega_count} / {total_points}) * {domain_area:.1f}")
    print(f"m(Omega) ≈ {omega_count / total_points:.4f} * {domain_area:.1f}")
    print(f"Estimated m(Omega) = {estimated_measure:.4f}")

if __name__ == '__main__':
    main()