import numpy as np

# The normalised density of the invariant measure is a linear function of the form rho(x) = a*x + b.
# For the map T(x) = 1/sqrt(x) mod 1, the density is given by rho(x) = 0.5 + x.

# The coefficient of x is 'a'.
a = 1.0

# The constant term is 'b'.
b = 0.5

print("The normalised density of the invariant measure is of the form: rho(x) = a*x + b")
print("The final equation for the density is:")
print(f"rho(x) = {a}*x + {b}")
print("The numbers in the final equation are:")
print(f"a = {a}")
print(f"b = {b}")
