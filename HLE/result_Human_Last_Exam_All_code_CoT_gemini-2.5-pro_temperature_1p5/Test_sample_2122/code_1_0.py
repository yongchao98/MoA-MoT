import numpy as np
from scipy.special import erfi

# We need to find T(t) at t_final, where t_final = arccosh(2)/2.
# This corresponds to x_final = cosh(2*t_final) = 2.
# The problem boils down to finding T(x=2).
# The final expression for T(x=2) is T(2) = -2 + 1/v(2),
# where v(2) is given by the expression below.

e = np.e
pi = np.pi

# Calculate the values needed for v(2)
# v(x) = exp(-x^2) * (e + (sqrt(pi)/2) * (erfi(x) - erfi(1)))
sqrt_pi_div_2 = np.sqrt(pi) / 2
erfi_of_1 = erfi(1)
erfi_of_2 = erfi(2)

# Calculate v(2)
v_of_2 = np.exp(-4) * (e + sqrt_pi_div_2 * (erfi_of_2 - erfi_of_1))

# Now calculate the final temperature T(2)
T_at_final_time = -2 + 1 / v_of_2

print(f"The equation for the temperature T at x=2 is: T(2) = -2 + 1/v(2)")
print(f"The equation for v(2) is: v(2) = exp(-4) * (e + (sqrt(pi)/2) * (erfi(2) - erfi(1)))")
print(f"e = {e}")
print(f"sqrt(pi)/2 = {sqrt_pi_div_2}")
print(f"erfi(1) = {erfi_of_1}")
print(f"erfi(2) = {erfi_of_2}")
print(f"exp(-4) = {np.exp(-4)}")
print(f"v(2) = {v_of_2}")
print(f"The temperature T at the specified time is T = {-2 + 1/v_of_2}")

<<<1.0863901615175924>>>