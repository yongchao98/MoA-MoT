import numpy as np
from scipy.special import erf

# The problem is to find the temperature T at a specific time t_f.
# The differential equation can be transformed and solved analytically,
# yielding a solution in terms of the error function (erf).

# The specific time is t_f = arccosh(2) / 2.
# We use the transformation x = cosh(2*t).
# At t_f, x = cosh(2 * arccosh(2) / 2) = cosh(arccosh(2)) = 2.
x = 2

# The analytical solution for the transformed equation with the initial condition T(t=0)=0 (which is T(x=1)=0) is:
# T(x) = -x + 1/v(x)
# where v(x) is given by:
# v(x) = exp(x^2 - 1) - exp(x^2) * (sqrt(pi)/2) * (erf(x) - erf(1))

# Calculate pi
pi_val = np.pi

# Calculate v(x=2)
# erf is the error function from scipy.special
v_at_2 = np.exp(2**2 - 1) - np.exp(2**2) * (np.sqrt(pi_val)/2) * (erf(2) - erf(1))

# Calculate T(x=2)
T_at_2 = -2 + 1/v_at_2

# The time t_f
t_f = np.arccosh(2) / 2

# Print the final result in the required format
print("The temperature T at time t = arccosh(2)/2 is:")
# The final equation is T(t_f) = T_at_2
# We print the numerical values for t_f and T_at_2
print(f"T({t_f}) = {T_at_2}")
