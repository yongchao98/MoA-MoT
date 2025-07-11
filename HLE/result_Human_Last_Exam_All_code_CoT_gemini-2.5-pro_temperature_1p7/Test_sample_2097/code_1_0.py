import math

# Based on the analysis, the minimum magnetization M_z(1) is achieved for n=1.
# The value is given by the expression A_1 = -32 / pi^2.

numerator = -32
pi_squared = math.pi**2
result = numerator / pi_squared

print(f"The minimum magnetization M_z(1) is found for n_min = 1.")
print(f"The calculation is: {numerator} / pi^2")
print(f"Which evaluates to: {numerator} / {pi_squared}")
print(f"The final numerical result is: {result}")