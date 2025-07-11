import numpy as np
from scipy.special import gamma

# Step 1: Define the initial value y(0) from the problem statement.
# y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
cbrt_3 = 3**(1/6)
gamma_2_3 = gamma(2/3)
y0 = 1 / (128 * cbrt_3 * gamma_2_3)

# Step 2: Calculate sec(pi/8).
# We use the half-angle identity for cosine: cos(x) = 2*cos^2(x/2) - 1
# Let x = pi/4, so x/2 = pi/8.
# cos(pi/4) = 2*cos^2(pi/8) - 1
# cos^2(pi/8) = (cos(pi/4) + 1) / 2
# cos(pi/4) = 1/sqrt(2)
cos_sq_pi_8 = (1/np.sqrt(2) + 1) / 2
# sec^2(pi/8) = 1 / cos^2(pi/8)
sec_sq_pi_8 = 2 / (1/np.sqrt(2) + 1)
# Take the square root to find sec(pi/8)
sec_pi_8 = np.sqrt(sec_sq_pi_8)

# Step 3: Calculate y(pi/4) using the assumed solution y(t) = y(0) * sec(t/2)
y_pi_4 = y0 * sec_pi_8

# As requested, printing the numbers in the final equation y(pi/4) = y(0) * sec(pi/8)
print(f"The initial value y(0) is: {y0}")
print(f"The value of sec(pi/8) is: {sec_pi_8}")
print(f"The final radius y(pi/4) is y(0) * sec(pi/8):")
print(y_pi_4)