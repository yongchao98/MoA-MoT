import numpy as np
from scipy.special import airy

# The constant K in the solution y(t) = K * (sqrt(3)Ai(sqrt(t)) + Bi(sqrt(t)))
K = 1 / 256.0

# The point at which to evaluate the solution
t = np.pi / 4
z = np.sqrt(t)

# Get the Airy function values Ai(z) and Bi(z)
# airy(z) returns (Ai, Aip, Bi, Bip)
ai_val, _, bi_val, _ = airy(z)

# The expression for the solution
sqrt3 = np.sqrt(3)
y_val = K * (sqrt3 * ai_val + bi_val)

# Output the components of the final equation
print(f"The solution is of the form y(t) = K * (sqrt(3) * Ai(sqrt(t)) + Bi(sqrt(t)))")
print(f"K = 1/256 = {K}")
print(f"We need to evaluate at t = pi/4, so z = sqrt(pi/4) = {z:.8f}")
print(f"Ai(z) = {ai_val:.8f}")
print(f"Bi(z) = {bi_val:.8f}")
print(f"The final equation is:")
print(f"y(pi/4) = {K} * (sqrt(3) * {ai_val:.8f} + {bi_val:.8f})")
print(f"y(pi/4) = {K} * ({sqrt3 * ai_val:.8f} + {bi_val:.8f})")
print(f"y(pi/4) = {K} * ({sqrt3 * ai_val + bi_val:.8f})")
print(f"y(pi/4) = {y_val:.8f}")
