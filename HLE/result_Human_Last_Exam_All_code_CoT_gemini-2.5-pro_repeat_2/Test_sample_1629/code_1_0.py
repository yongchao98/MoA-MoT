import numpy as np
from scipy.integrate import solve_ivp

def ode_system(t, y):
    """Defines the system of differential equations."""
    a, b = y
    a_prime = -b * a
    b_prime = -b**2 / 2 - a**2 + 6 * a - 6
    return [a_prime, b_prime]

def event_blowup(t, y):
    """Event function to detect when solution leaves a large box."""
    # Stop integration if |a| or |b| exceeds a large threshold.
    # This indicates a potential blow-up.
    return max(abs(y[0]), abs(y[1])) - 10000.0
event_blowup.terminal = True
event_blowup.direction = 1

# Define the grid of initial conditions
a_min, a_max = -1.0, 1.0
b_min, b_max = 2.0, 3.0
na, nb = 101, 51  # Number of points in each dimension
a_points = np.linspace(a_min, a_max, na)
b_points = np.linspace(b_min, b_max, nb)

blowup_count = 0
total_points = na * nb

# Iterate over the grid of initial conditions
for a0 in a_points:
    for b0 in b_points:
        y0 = [a0, b0]
        # Solve the ODE
        sol = solve_ivp(ode_system, [0, 20], y0, events=event_blowup, method='RK45')

        # Check if the integration was terminated by the event
        if sol.t_events[0].size > 0:
            a_final, b_final = sol.y_events[0][0]
            # Check if it's the correct type of blow-up (a -> +inf, b -> -inf)
            if a_final > 0 and b_final < 0:
                blowup_count += 1

# Calculate the measure of the set Omega
area_D = (a_max - a_min) * (b_max - b_min)
fraction_blowup = blowup_count / total_points
m_omega_estimated = fraction_blowup * area_D

# Output the calculation
print(f"The system was simulated for {total_points} initial conditions.")
print(f"Number of points leading to the specified blow-up: {blowup_count}")
print("The measure of the set Omega is estimated by the equation:")
print(f"m(Omega) = (blow-up count / total points) * Area(D)")
print(f"m(Omega) = ({blowup_count} / {total_points}) * {area_D:.1f}")
print(f"m(Omega) ~= {fraction_blowup:.4f} * {area_D:.1f} = {m_omega_estimated:.4f}")
print(f"\nThe estimated measure m(Omega) is approximately: {round(m_omega_estimated)}")
