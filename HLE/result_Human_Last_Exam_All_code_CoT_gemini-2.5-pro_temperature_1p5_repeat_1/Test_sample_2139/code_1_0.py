import numpy as np
from scipy.special import gamma

# The solution to the differential equation, as verified numerically,
# is y(t) = y(0) * cos(t^2).
# We need to calculate y(pi/4).

# 1. Calculate the initial value y(0).
c1 = 128.0
c2 = 3**(1/6)
c3 = gamma(2/3)
y0_inv = c1 * c2 * c3
y0 = 1 / y0_inv

# 2. Calculate the value of the function f(t) = cos(t^2) at t = pi/4.
t_final = np.pi / 4
cos_val = np.cos(t_final**2)

# 3. The final result is the product of y(0) and f(pi/4).
y_final = y0 * cos_val

# As requested, output the numbers in the final equation.
# The final equation is y(pi/4) = y(0) * cos((pi/4)^2)
print(f"The initial value y(0) is calculated as 1 / (c1 * c2 * c3)")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"So, y(0) = {y0}")
print(f"The value of cos((pi/4)^2) is calculated.")
print(f"pi/4 = {t_final}")
print(f"(pi/4)^2 = {t_final**2}")
print(f"cos((pi/4)^2) = {cos_val}")
print(f"The final result y(pi/4) = y(0) * cos((pi/4)^2)")
print(f"y(pi/4) = {y0} * {cos_val}")
print(f"y(pi/4) = {y_final}")
