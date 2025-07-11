import math

# The value of the integral is given by the formula b*f(b) - a*f(a).
a = 0
b = 2

# We calculated f(b) = f(2) = 2**(1/16).
f_b = 2**(1/16)

# We calculated f(a) = f(0) = 0.
f_a = 0

# The final result is b*f(b) - a*f(a).
result = b * f_b - a * f_a

print("The integral is evaluated using the inverse function integration property.")
print("Value = b*f(b) - a*f(a)")
print(f"a = {a}")
print(f"b = {b}")
print(f"f(a) = {f_a}")
print(f"f(b) = 2^(1/16)")
print("Result = 2 * 2^(1/16) - 0 * 0")
print("Result = 2^(1) * 2^(1/16)")
print("Result = 2^(1 + 1/16)")
print("Result = 2^(17/16)")
print(f"The numerical value is: {result}")
