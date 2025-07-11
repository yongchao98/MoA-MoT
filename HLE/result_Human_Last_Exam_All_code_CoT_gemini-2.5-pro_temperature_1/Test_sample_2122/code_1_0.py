import numpy as np

# Step 1: Define the time t_f
# t_f = arccosh(2) / 2
# Let u = 2 * t_f = arccosh(2)
# By definition, cosh(u) = 2

# Step 2: Calculate sinh(u) using the identity cosh^2(u) - sinh^2(u) = 1
# sinh^2(u) = cosh^2(u) - 1 = 2^2 - 1 = 3
# sinh(u) = sqrt(3)
sinh_2tf = np.sqrt(3)

# Step 3: The solution to the differential equation at the specific point is T(t_f) = tanh(sinh(2*t_f))
# So, T(t_f) = tanh(sqrt(3))
result = np.tanh(sinh_2tf)

# Step 4: Print the final equation and the result.
# The problem asks to output each number in the final equation.
# The final equation is T = tanh(sqrt(3))
number_in_equation = 3
print(f"The final equation is T = tanh(sqrt({number_in_equation}))")
print(f"The temperature at the specified time is: {result}")
