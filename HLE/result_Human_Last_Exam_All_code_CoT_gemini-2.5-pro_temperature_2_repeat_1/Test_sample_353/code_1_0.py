import numpy as np
from scipy.optimize import minimize_scalar

# 1. Define the stability boundary curve z(theta) for BDF4.
# z(theta) = 25/12 - 4*exp(-i*theta) + 3*exp(-2i*theta) - 4/3*exp(-3i*theta) + 1/4*exp(-4i*theta)
def get_z(theta):
    """Calculates the complex value on the stability boundary for a given theta."""
    return (25/12 
            - 4 * np.exp(-1j * theta) 
            + 3 * np.exp(-2j * theta) 
            - (4/3) * np.exp(-3j * theta) 
            + (1/4) * np.exp(-4j * theta))

# 2. Define the function to maximize, which is the angle of z(theta).
# We minimize its negative to find the maximum.
def angle_to_maximize(theta):
    """Returns the negative angle of z(theta). For use with a minimizer."""
    z = get_z(theta)
    return -np.arctan2(z.imag, z.real)

# 3. Use a numerical optimizer to find the theta that maximizes the angle.
# We search in the interval [0, pi] due to symmetry.
opt_result = minimize_scalar(angle_to_maximize, bounds=(0, np.pi), method='bounded')
theta_max = opt_result.x
phi_max = -angle_to_maximize(theta_max)

# 4. Calculate the A(alpha)-stability angle alpha.
alpha = np.pi - phi_max

# 5. Determine the integer N for the exact expression arctan(sqrt(N)).
# We expect tan(alpha)^2 to be an integer.
val = np.tan(alpha)**2
int_val = int(round(val))

# 6. Print the results step-by-step.
print(f"Numerical search for the angle of tangency:")
print(f"  - The angle of the curve z(theta) is maximized at theta = {np.rad2deg(theta_max):.4f} degrees.")
print(f"  - The maximum angle of the curve, phi_max, is {np.rad2deg(phi_max):.4f} degrees.")
print(f"\nThe A(alpha)-stability angle is given by the formula alpha = pi - phi_max.")
print(f"  - alpha = 180 - {np.rad2deg(phi_max):.4f} = {np.rad2deg(alpha):.4f} degrees.")
print(f"\nTo find the exact expression, we compute tan(alpha)^2:")
print(f"  - tan({np.rad2deg(alpha):.4f})^2 = {val:.4f}, which is very close to the integer {int_val}.")
print(f"\nThis suggests the exact value for the angle is alpha = arctan(sqrt(value)).")
print(f"Thus, the final equation for alpha is:")
print(f"alpha = arctan(sqrt({int_val}))")