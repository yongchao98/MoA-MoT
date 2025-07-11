import math

# Coefficients from the PDE:
# \frac{\partial u}{\partial t} + alpha * u * \frac{\partial u}{\partial x} - beta * \frac{\partial^2 u}{\partial x^2} = 0
alpha = 4
beta = 3

# Given condition for the spatial gradient at the stationary point x_0 at time tau
u_x = -1.0

# The value of the integral can be found using a specific identity derived for
# this type of solution. The identity relates the integral to the coefficients
# of the PDE and the local gradient.
# I = (3 / 20) * alpha^2 * beta * |u_x|^3
# Let's break down the calculation.
term1 = 3
term2 = 20
term3 = alpha**2
term4 = beta
term5 = abs(u_x)**3

# Calculate the final result
integral_value = (term1 / term2) * term3 * term4 * term5

# Output the equation with the final value.
print(f"The integral is calculated by the formula (3/20) * alpha^2 * beta * |u_x|^3")
print(f"Plugging in the values alpha={alpha}, beta={beta}, u_x={u_x}:")
print(f"Result = ({term1} / {term2}) * {alpha}^2 * {beta} * |{u_x}|^3")
print(f"Result = ({term1} / {term2}) * {term3} * {term4} * {term5}")
print(f"Result = {integral_value}")
