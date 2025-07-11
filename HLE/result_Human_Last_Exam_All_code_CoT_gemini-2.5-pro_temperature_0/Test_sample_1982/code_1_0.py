import math

# Given values
c1 = 1e4

# Expressions for l1 and alpha1
# l1 = (1 + 10^5)^5
# alpha1 = (1 + 10^5)^6 * (1 - 10^5 + 10^10)

# As derived in the explanation, the value of x11 can be found by
# x11 = alpha1 / l1
# x11 = (1 + 1e5) * (1 - 1e5 + 1e10)
# This is of the form (a+b)(a^2 - ab + b^2) = a^3 + b^3
# with a=1 and b=1e5.
a = 1.0
b = 1e5
x11 = a**3 + b**3

# The key equation derived from the matrix problem is x11 = 1 + c1 * u1
# We can solve for u1: u1 = (x11 - 1) / c1
u1 = (x11 - 1) / c1

# Print the final equation with the calculated numbers
# The equation is x11 = 1 + c1 * u1
print("Derived equation with numerical values:")
print(f"{x11} = 1 + {c1} * {u1}")

# Print the final answer for u1
print("\nThe value of the control u1 is:")
print(u1)