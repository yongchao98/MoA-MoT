import math

# Based on the derivation using the Local Limit Theorem, the limit is
# (Number of points) * (Area of lattice cell) * (Value of PDF at origin).
# Number of points = 3
# Area of lattice cell = sqrt(3)
# PDF of the approximating Normal N(0, n/2 * I) at the origin is 1/(pi*n).
# So, P(n) is approximately 3 * sqrt(3) / (pi * n).
# The limit of n * P(n) as n tends to infinity is 3 * sqrt(3) / pi.

# The numbers in the final equation for the limit are 3, sqrt(3), and pi.
num1 = 3
num2_val = math.sqrt(3)
num2_str = "sqrt(3)"
den_val = math.pi
den_str = "pi"

result = (num1 * num2_val) / den_val

print(f"The final expression for the limit is: ({num1} * {num2_str}) / {den_str}")
print(f"The number {num1} is an integer: {num1}")
print(f"The number {num2_str} is the square root of 3: {num2_val}")
print(f"The number {den_str} is the mathematical constant pi: {den_val}")
print(f"The value of the limit is {result}")