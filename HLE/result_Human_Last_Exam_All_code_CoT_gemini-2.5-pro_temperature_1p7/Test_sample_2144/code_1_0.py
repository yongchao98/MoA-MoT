import math

# The problem is to find x0 for y(x0)=-3 given the ODE and y(0)=-1.
# The analytical solution for the particle's trajectory is:
# 3 * y**(1/3) = (1/2)**(1/3) * x - 3

# Define the target y value
y_target = -3

# Define the numbers in the final equation we need to solve for x0
# Equation is: 3 * (-3)**(1/3) = (1/2)**(1/3) * x0 - 3
coeff_y = 3.0
base_y = -3.0
power_y = 1.0/3.0
base_x_denom = 2.0
const_term = -3.0

print(f"The equation for the particle's trajectory is: {coeff_y} * y^({power_y}) = (1/{base_x_denom})^({power_y}) * x + {const_term}")
print(f"Solving for x0 when y = {base_y}:")
print(f"{coeff_y} * ({base_y})^({power_y}) = (1/{base_x_denom})^({power_y}) * x0 + {const_term}")

# To compute the real cube root of a negative number `n`, we calculate `-abs(n)**(1/3)`.
lhs = coeff_y * (-math.pow(abs(base_y), power_y))

# The coefficient of x0 on the right hand side
rhs_coeff = math.pow(1.0/base_x_denom, power_y)

# Solve for x0: lhs = rhs_coeff * x0 + const_term
x0 = (lhs - const_term) / rhs_coeff

print(f"\nThe calculated position x0 is: {x0}")