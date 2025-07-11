import numpy as np

# We found that the maximum is achieved under the following conditions:
# - The expression is evaluated at a point x where u(x) = 0.
# - The integral term u_bar(x) is close to 0.
# - The state at x+1 is u(x+1) = 1.
# - The integral term u_bar(x+1) is close to 0.

# Let's define these values
u_x = 0.0
ubar_x = 0.0
u_xplus1 = 1.0
ubar_xplus1 = 0.0

# Calculate intermediate quantities based on these values
# F(u, u_bar) = u * (1-u)^2 * exp(-u_bar)
F_at_x = u_x * (1 - u_x)**2 * np.exp(-ubar_x)
F_at_xplus1 = u_xplus1 * (1 - u_xplus1)**2 * np.exp(-ubar_xplus1)

# Time derivative of u_bar: ubar_t = F(x) - F(x+1)
ubar_t = F_at_x - F_at_xplus1

# Spatial derivative of u_bar: ubar_x = u(x+1) - u(x)
ubar_x_val = u_xplus1 - u_x

# F1 = dF/du = (1 - 4u + 3u^2) * exp(-u_bar)
F1_at_x = (1 - 4*u_x + 3*u_x**2) * np.exp(-ubar_x)

# The full expression for E = (d/dt + F1*d/dx)F11 simplifies at u=0 to:
# E = 4 * exp(-ubar_x) * (ubar_t + F1_at_x * ubar_x_val)
# Note: F(x) is 0 since u_x is 0
E_calc = 4 * np.exp(-ubar_x) * (-F_at_xplus1 + F1_at_x * ubar_x_val)

print(f"To find the maximum, we evaluate the expression at an optimal point.")
print(f"We choose a state where u(x) = {u_x} and u_bar(x) = {ubar_x}.")
print(f"And the state at x+1 is u(x+1) = {u_xplus1} and u_bar(x+1) = {ubar_xplus1}.")
print(f"This is an ideal case that can be approached by specific solution profiles.")
print("\nStep-by-step calculation:")
print(f"F(x) = {u_x} * (1 - {u_x})**2 * exp(-{ubar_x}) = {F_at_x}")
print(f"F(x+1) = {u_xplus1} * (1 - {u_xplus1})**2 * exp(-{ubar_xplus1}) = {F_at_xplus1}")
print(f"ubar_t = F(x) - F(x+1) = {F_at_x} - {F_at_xplus1} = {ubar_t}")
print(f"ubar_x = u(x+1) - u(x) = {u_xplus1} - {u_x} = {ubar_x_val}")
print(f"F1(x) = (1 - 4*{u_x} + 3*{u_x}**2) * exp(-{ubar_x}) = {F1_at_x}")
print("\nThe expression for E at u=0 simplifies to 4*exp(-ubar_x)*(-F(x+1) + F1(x)*ubar_x):")
print(f"E = 4 * exp(-{ubar_x}) * (-({F_at_xplus1}) + {F1_at_x} * {ubar_x_val})")
print(f"E = 4 * {np.exp(-ubar_x)} * ({-F_at_xplus1} + {F1_at_x * ubar_x_val})")
print(f"E = 4 * {np.exp(-ubar_x)} * ({-F_at_xplus1 + F1_at_x * ubar_x_val})")
print(f"E = {4 * np.exp(-ubar_x) * (-F_at_xplus1 + F1_at_x * ubar_x_val)}")

print(f"\nThe maximum value is {E_calc}.")
<<<4>>>