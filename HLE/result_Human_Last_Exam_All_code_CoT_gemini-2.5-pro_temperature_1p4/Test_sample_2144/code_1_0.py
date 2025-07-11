import numpy as np

# Step 1: State the parametric solution for the trajectory.
# The ODE can be solved using a parametric method. The particle's trajectory
# is found to be described by the following equations involving a parameter 'c':
#   1) y = c**2 * (x - c)
#   2) x**2 - 4*c*x + 3*c**2 = K  (where K is a constant)

# Step 2: Use the initial condition y(0) = -1 to find K and the initial value of c (c0).
print("Step 1: Applying the initial condition y(0) = -1.")
# At x = 0, y = -1. From equation (1):
# -1 = c0**2 * (0 - c0)  =>  -1 = -c0**3  =>  c0**3 = 1
# The only real solution is c0 = 1.
c0 = 1.0
print(f"The parameter 'c' at x=0 is c0 = {c0}.")

# Now, use c0=1 and x=0 in equation (2) to find K:
# 0**2 - 4*(1)*0 + 3*(1)**2 = K  =>  K = 3
K = 3.0
print(f"The constant of integration K for this trajectory is {K}.")
print("-" * 30)


# Step 3: Find the value of parameter 'c' at the target point where y = -3.
print("Step 2: Finding the parameter 'c' where y = -3.")
# We have a system of two equations for x and c:
#   1) x**2 - 4*c*x + 3*c**2 = 3
#   2) -3 = c**2 * (x - c)

# From equation (2), we can express x as: x = c - 3/c**2.
# Substituting this into equation (1) yields a polynomial equation for c:
# c**4 - 2*c**3 - 3 = 0
print("This leads to solving the polynomial equation: c^4 - 2*c^3 - 3 = 0")

# Coefficients of the polynomial: 1*c^4 - 2*c^3 + 0*c^2 + 0*c - 3 = 0
coeffs = [1, -2, 0, 0, -3]
roots = np.roots(coeffs)

# The roots can be complex. We are interested in the real roots.
real_roots = roots[np.isreal(roots)].real
print(f"The real roots for 'c' are: {real_roots}")

# The trajectory starts at x=0 with c=1. The value of c changes as the particle moves.
# The root c = -1.0 is the correct one for the continuous path from the initial state.
c_target = -1.0
print(f"Selecting the relevant root for the trajectory: c = {c_target}.")
print("-" * 30)


# Step 4: Calculate the position x0 using the found value of c.
print("Step 3: Calculating the position x0.")
# Use the relation x = c - 3/c**2 derived earlier.
print(f"The formula for x0 is: x0 = c - 3 / c^2")
print(f"Substituting c = {c_target}:")
# The final equation with numbers
print(f"x0 = {c_target} - 3 / ({c_target})^2")
c_target_sq = c_target**2
print(f"x0 = {c_target} - 3 / {c_target_sq}")
div_result = 3 / c_target_sq
print(f"x0 = {c_target} - {div_result}")
x0 = c_target - div_result
print(f"x0 = {x0}")

# Final Answer
# print(f"\n<<<The position x0 where y=-3 is {x0}.>>>")