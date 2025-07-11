import numpy as np

# Step 1: Under the hypothesis that det(N) = 0, the characteristic function is 1.
# The integral for phi(a) simplifies significantly. We need to evaluate phi(7).
a = 7

# Step 2: The simplified integral is phi(a) = integral from 0 to infinity of -2 * sin(a*t)/t dt.
# This is -2 times the Dirichlet integral.
dirichlet_integral_val = np.pi / 2

# Step 3: Calculate the value of phi(7).
# For a=7 > 0, the result is -2 * sgn(7) * (pi/2).
sgn_a = np.sign(a)
phi_7 = -2 * sgn_a * dirichlet_integral_val

# Step 4: Output the equation and the final value.
print(f"The equation simplifies to: phi(7) = integral from 0 to inf of (-2 * sin(7*t) / t) dt")
print(f"This evaluates to: -2 * sgn(7) * (pi / 2)")
print(f"                                = {-2} * {int(sgn_a)} * ({np.pi} / {2})")
print(f"                                = {phi_7}")
print("\nThe exact value of phi(7) is -pi.")
print(f"Final Value: {phi_7}")
