import scipy.integrate as integrate
import numpy as np

# Define the value of z
z = 0.2

# Integrand for Region 1 where y >= 2x
def integrand_r1(y, x):
    # This region corresponds to p = y
    # The term is 6 * p * (1-p) / (y-x)
    return 6 * y * (1 - y) / (y - x)

# Integrand for Region 2 where y < 2x
def integrand_r2(y, x):
    # This region corresponds to p = 2 * (y-x)
    # The term is 6 * p * (1-p) / (y-x) which simplifies
    return 12 * (1 - 2 * (y - x))

# --- Numerical Integration ---

# The integration domain for Region 1 is split into two parts
# Part 1a: x from 0 to z/2 (i.e., 0.1), y from z to 1
# In this part, y > z >= 2x, so we are always in Region 1
integral_1a, _ = integrate.dblquad(integrand_r1, 0, z / 2, lambda x: z, lambda x: 1)

# Part 1b: x from z/2 to z (i.e., 0.1 to 0.2), y from 2x to 1
integral_1b, _ = integrate.dblquad(integrand_r1, z / 2, z, lambda x: 2 * x, lambda x: 1)

# The integration domain for Region 2
# x from z/2 to z (i.e., 0.1 to 0.2), y from z to 2x
integral_2, _ = integrate.dblquad(integrand_r2, z / 2, z, lambda x: z, lambda x: 2 * x)

# The total integral is the sum of these parts
total_integral_one_side = integral_1a + integral_1b + integral_2

# The PDF f_Z(z) is twice this value due to symmetry
f_z_value = 2 * total_integral_one_side

print(f"The calculation is f({z}) = 2 * (I_1a + I_1b + I_2)")
print(f"where I_1a = {integral_1a}")
print(f"      I_1b = {integral_1b}")
print(f"      I_2  = {integral_2}")
print(f"f({z}) = 2 * ({integral_1a} + {integral_1b} + {integral_2})")
print(f"f({z}) = 2 * ({total_integral_one_side})")
print(f"The final value is f({z}) = {f_z_value}")
<<<3.6>>>