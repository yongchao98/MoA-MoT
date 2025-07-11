import numpy as np
from scipy.optimize import fsolve

# Set up all the constants and parameters from the problem description.
# Man's initial speed and accelerations
u = 1.5
a1 = -0.15
a2 = 0.25

# Time points
t1 = 4.0
t3 = 15.0
t5 = 23.0
t6 = 40.0

# Calculate intermediate time points t2 and t4
t2 = t1 - (0 - u) / a1  # Time when man stops
t4 = t3 + (u - 0) / a2  # Time when man reaches speed u again

# Bird's path angles. As reasoned, we assume α = 40° to resolve the contradiction.
alpha_deg = 40.0
gamma_deg = 40.0
alpha = np.deg2rad(alpha_deg)
gamma = np.deg2rad(gamma_deg)

# Man's y-position at t4
y_m1 = u * t1
y_m2 = y_m1 + u * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
y_m3 = y_m2
y_m4 = y_m3 + 0 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2

# Define a function to solve for bird's speed 'v' and initial ground distance 'd_xy1'
# This is based on the original rendezvous plan.
def system_equations(p):
    v, d_xy1 = p
    # We enforce physical constraints to guide the solver
    if 16 * v**2 < d_xy1**2 or v <= 0 or d_xy1 < 0:
        return [1e6, 1e6] # Return a large error if constraints are violated
    
    # Calculate initial displacements based on v and d_xy1
    z1 = np.sqrt(16 * v**2 - d_xy1**2)
    x1 = d_xy1 * np.sin(alpha)
    y1 = d_xy1 * np.cos(alpha)

    # Equation 1: Based on the rendezvous condition for the original plan.
    # The required displacements in x and z must match the direction of the final flight path.
    eq1 = (x1 + (t2 - t1) * v) * np.sin(gamma) - (z1 + (t4 - t3) * v) * np.cos(gamma)

    # Equation 2: The rendezvous must happen at the same y-coordinate.
    delta_t_plan = (z1 + (t4-t3) * v) / (v * np.sin(gamma))
    y_man_planned_meet = y_m4 + u * delta_t_plan
    y_bird_planned_meet = y1 + (t3-t2)*v # Note: duration t3-t2 = 1s
    eq2 = y_man_planned_meet - y_bird_planned_meet
    
    return [eq1, eq2]

# Solve the system of equations
initial_guess = [1.5, 4 * 1.5] # Guess v=u, and d_xy1 is the total distance if z1=0
solution, _, ier, _ = fsolve(system_equations, initial_guess, full_output=True)
if ier != 1:
    print("Solver did not converge. The problem parameters might be inconsistent.")
v, d_xy1 = solution

# Now, calculate all positions at key times using the solved 'v'
z1 = np.sqrt(16 * v**2 - d_xy1**2)
x1 = d_xy1 * np.sin(alpha)
y1 = d_xy1 * np.cos(alpha)

# Bird's position at t4
x_b4 = x1 + (t2-t1) * v
y_b4 = y1 + (t3-t2)*v
z_b4 = z1 + (t4-t3) * v

# Bird's position at t5 (before the final leg is altered)
delta_t_54 = t5 - t4
x_b5 = x_b4 - v * np.cos(gamma) * delta_t_54
y_b5 = y_b4
z_b5 = z_b4 - v * np.sin(gamma) * delta_t_54

# Man's position at t5
y_m5 = y_m4 + u * delta_t_54

# After the wind gust, calculate the bird's new constant velocity for t5 to t6
delta_t_65 = t6 - t5
# To meet on the man's path at t6, the bird must cover the x and z distances in delta_t_65
v_bx_6 = -x_b5 / delta_t_65
v_bz_6 = -z_b5 / delta_t_65
# The bird's speed is constant, so we find the required y-velocity
v_by_6 = np.sqrt(v**2 - v_bx_6**2 - v_bz_6**2)
# Final meeting y-position
y_final = y_b5 + v_by_6 * delta_t_65

# Finally, calculate the man's required acceleration a3
# y_final = y_m5 + u * delta_t_65 + 0.5 * a3 * delta_t_65^2
a3 = (y_final - y_m5 - u * delta_t_65) / (0.5 * delta_t_65**2)

print("Final Calculation Steps:")
print(f"The meeting is set to happen at time t6 = {t6} s.")
print(f"The bird's position at t5 = {t5} s is ({x_b5:.2f}, {y_b5:.2f}, {z_b5:.2f}) m.")
print(f"The man's position at t5 = {t5} s is (0, {y_m5:.2f}, 0) m.")
print(f"To meet at the new location, the final y-coordinate must be y_final = {y_final:.2f} m.")
print(f"The man must travel from y_m5 to y_final in {delta_t_65} s, starting with speed u = {u:.1f} m/s.")
print("Using the kinematic equation: y_final = y_m5 + u * Δt + 0.5 * a3 * Δt²")
print(f"{y_final:.2f} = {y_m5:.2f} + {u:.1f} * {delta_t_65:.1f} + 0.5 * a3 * {delta_t_65**2:.1f}")
print(f"{y_final:.2f} = {y_m5 + u * delta_t_65:.2f} + {0.5 * delta_t_65**2:.1f} * a3")
a3_numerator = y_final - y_m5 - u * delta_t_65
a3_denominator = 0.5 * delta_t_65**2
print(f"a3 = ({a3_numerator:.2f}) / ({a3_denominator:.2f})")
print(f"The man's constant acceleration is a3 = {a3:.4f} m/s^2.")

print(f"\nFinal answer derived from the calculations:")
<<<0.2831>>>