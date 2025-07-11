import numpy as np
import sympy

# Part 1: Calculation of the sum
# The problem defines l(n,p) as the injectivity radius of the Stiefel manifold M(n,p).
# For n > p, this value is pi. The primes p_(k) used in the sum satisfy this.
# The sum is Sum_{i=1 to 10} Sum_{j=1 to 10} pi = 10 * 10 * pi = 100 * pi.
sum_val = 100 * np.pi

# Part 2: Calculation of the integral
# The integral simplifies to Integral_g + Integral(x*e^(-x) dx)
# Integral(x*e^(-x) dx) from 0 to inf is 1.

# For Integral_g, we first calculate the dimensions d1 and d2.
# dim(M(n, p)) = n*p - p*(p+1)/2

# Get prime numbers using sympy
n1_idx = 8231
p1_idx = 781
n1 = sympy.prime(n1_idx)
p1 = sympy.prime(p1_idx)

n2_idx = 10231
p2_idx = 2321
n2 = sympy.prime(n2_idx)
p2 = sympy.prime(p2_idx)

# Using integer arithmetic for precision
d1 = n1 * p1 - (p1 * (p1 + 1)) // 2
d2 = n2 * p2 - (p2 * (p2 + 1)) // 2

# The exponents d1 and d2 are extremely large.
# As explained in the plan, the integral term involving these exponents, I_g, tends to 0.
# The term 1/(1+x^(2d)) becomes a step function for large d.
# The difference of two such functions for d1 and d2 makes the integral evaluate to 0.
integral_g = 0
integral_val = integral_g + 1

# Final calculation
result = sum_val * integral_val

print(f"Step 1: Calculate the dimensions d1 and d2.")
print(f"p_({n1_idx}) = {n1}, p_({p1_idx}) = {p1}")
print(f"d1 = dim(M({n1}, {p1})) = {d1}")
print(f"p_({n2_idx}) = {n2}, p_({p2_idx}) = {p2}")
print(f"d2 = dim(M({n2}, {p2})) = {d2}")
print("-" * 20)
print("Step 2: Calculate the two main factors of the expression.")
print(f"The first factor (the sum) is 100 * pi.")
print(f"Value of the sum = {sum_val}")
print("The second factor (the integral) evaluates to 1, as the complex part of the integral is assumed to be 0 due to the large dimensions.")
print(f"Value of the integral = {integral_val}")
print("-" * 20)
print("Step 3: Calculate the final result.")
print(f"Final Result = ({sum_val}) * ({integral_val})")
print(f"The final value is: {result}")
