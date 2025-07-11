import math
import cmath

# Based on the analysis, a sign in the PDE is assumed to be corrected for a consistent solution.
# u_t - 6*u*u_x + u_xxx - 5*u_xx = 0.
# This leads to a traveling wave solution u(x,t) = F(x-c*t) with wave speed c=6.

# 1. The function to be differentiated is f(x,t) = g(x-6t).
#    The fractional derivative operator D_t^(1/2) D_x^(1/2) simplifies to
#    multiplication by (-c)^(1/2) = (-6)^(1/2) = i*sqrt(6), and one derivative d/dz.
c = 6
operator_factor = cmath.sqrt(-c)
print(f"The wave speed is c = {c}")
print(f"The fractional derivative operator simplifies to a factor of {operator_factor:.6f}")

# 2. The transformed function is g(z) = (1+exp(-z))^2 / (2*sqrt(6)).
#    Its derivative is g'(z) = -(exp(-z) + exp(-2*z)) / sqrt(6).
# 3. We evaluate the final expression at x=6t, which makes the argument z = x-6t = 0.
z = 0
print(f"The expression is evaluated at z = {z}")

# 4. Calculate g'(z) at z=0.
g_prime_at_0 = -(math.exp(-z) + math.exp(-2*z)) / math.sqrt(c)
print(f"The value of the derivative g'(z) at z=0 is: {g_prime_at_0:.6f}")

# 5. Compute the final result: operator_factor * g'(0)
final_result = operator_factor * g_prime_at_0
print("\nThe final computed quantity is:")
print(final_result)