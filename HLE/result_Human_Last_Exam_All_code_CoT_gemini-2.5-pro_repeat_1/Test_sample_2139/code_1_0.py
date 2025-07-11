import numpy as np
from scipy.special import gamma

# The solution to the differential equation that satisfies the initial conditions
# is y(t) = y(0) * cos(t^2) * sqrt(cos(t)).
# We need to find the value of this function at t = pi/4.
# The final equation to compute is:
# y(pi/4) = y(0) * cos((pi/4)^2) * sqrt(cos(pi/4))

# Step 1: Calculate the initial value y(0) from the given expression.
# y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
gamma_val = gamma(2/3)
term_3_pow_1_6 = 3**(1/6)
y0 = (128 * term_3_pow_1_6 * gamma_val)**(-1)

# Step 2: Define the target time t and calculate the remaining terms of the solution.
t_final = np.pi / 4
cos_t_squared_val = np.cos(t_final**2)
sqrt_cos_t_val = np.sqrt(np.cos(t_final))

# Step 3: Calculate the final result by multiplying the components.
y_at_pi_over_4 = y0 * cos_t_squared_val * sqrt_cos_t_val

# As requested, output each number in the final equation.
print("The calculation is based on the equation: y(pi/4) = y(0) * cos((pi/4)^2) * sqrt(cos(pi/4))")
print(f"Value of y(0): {y0}")
print(f"Value of cos((pi/4)^2): {cos_t_squared_val}")
print(f"Value of sqrt(cos(pi/4)): {sqrt_cos_t_val}")
print(f"\nThe final calculated value y(pi/4) is the product of these three numbers:")
print(f"y(pi/4) = {y0} * {cos_t_squared_val} * {sqrt_cos_t_val}")
print(f"\nResult: {y_at_pi_over_4}")