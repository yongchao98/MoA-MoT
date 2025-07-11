import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

def solve_and_estimate_area():
    """
    Solves the given system of ODEs for a grid of initial conditions
    and estimates the area of the set Omega where solutions blow up.
    """
    # Define the system of ordinary differential equations
    def ode_system(t, y):
        a, b = y
        dadt = -b * a
        dbdt = -b**2 / 2 - a**2 + 6 * (a - 1)
        return [dadt, dbdt]

    # Define the domain for initial conditions
    a0_range = [-1, 1]
    b0_range = [2, 3]
    domain_area = (a0_range[1] - a0_range[0]) * (b0_range[1] - b0_range[0])

    # Set up a grid of initial points.
    # Using an even number for n_a ensures the grid is symmetric around a=0
    # and doesn't contain a=0 itself, leading to a more accurate estimate.
    n_a = 40
    n_b = 20
    a0_vals = np.linspace(a0_range[0], a0_range[1], n_a)
    b0_vals = np.linspace(b0_range[0], b0_range[1], n_b)
    total_count = n_a * n_b

    # Set up simulation parameters
    t_max = 20.0  # Maximum integration time
    blowup_threshold = 1000.0  # Threshold for 'a' to be considered as blowing up

    # Define an event to detect when a(t) becomes very large and positive
    def blowup_event(t, y):
        return y[0] - blowup_threshold
    blowup_event.terminal = True  # Stop integration when this event occurs
    blowup_event.direction = 1    # Trigger event only when a(t) is increasing

    # Counter for initial conditions in Omega
    omega_count = 0

    # Iterate over the grid of initial conditions
    for a0 in a0_vals:
        for b0 in b0_vals:
            initial_conditions = [a0, b0]
            
            # Solve the ODE for the current initial condition
            sol = solve_ivp(
                ode_system,
                [0, t_max],
                initial_conditions,
                events=blowup_event,
                dense_output=True
            )
            
            # Check if the blow-up event was triggered
            if sol.t_events[0].size > 0:
                omega_count += 1

    # Estimate the area of Omega
    estimated_area = (omega_count / total_count) * domain_area

    # Print the results, showing the numbers used in the final calculation
    print(f"Numerical Estimation of m(Omega):")
    print("-" * 35)
    print(f"Initial conditions domain R: [a,b] in [-1,1] x [2,3]")
    print(f"Total area of R: {domain_area}")
    print(f"Grid points sampled: {total_count} ({n_a} for 'a' and {n_b} for 'b')")
    print(f"Points in Omega (a -> +inf): {omega_count}")
    print(f"The equation for the area estimate is: m(Omega) = (omega_count / total_count) * domain_area")
    print(f"Plugging in the numbers: m(Omega) = ({omega_count} / {total_count}) * {domain_area}")
    print(f"Estimated size of the set Omega: {estimated_area:.4f}")

# Run the simulation and print the output
solve_and_estimate_area()