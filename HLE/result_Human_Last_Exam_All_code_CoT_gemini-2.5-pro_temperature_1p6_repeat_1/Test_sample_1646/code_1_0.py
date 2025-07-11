import sympy
from sympy import pi, gamma, Symbol, pprint

# Define the symbol for the complex variable z
z = Symbol('z')

# Based on the derivation, the function f(z) has a constant numerator
# and the Gamma function in the denominator.

# The constant in the numerator is the square root of pi.
constant_pi = sympy.sqrt(pi)

# The numerator is sqrt(pi).
# According to the instructions, we output each number in the final equation.
# The number involved is pi.
print(f"The number used in the numerator is pi: {float(pi)}")
print("The numerator of the function is the square root of pi:")
pprint(constant_pi)
print("-" * 20)

# The denominator is the Gamma function of z, denoted as Gamma(z) or gamma(z).
denominator_func = gamma(z)
print("The denominator of the function is Gamma(z):")
pprint(denominator_func)
print("-" * 20)

# The explicit form of f(z) is the ratio of the two.
f_z = constant_pi / denominator_func

# Print the final explicit form of f(z).
print("The final explicit form of f(z) is:")
pprint(f_z)
