import numpy as np

# The problem has been reduced to solving a system of algebraic equations.
# The core of this is finding the roots of the polynomial: z^4 - 2*z^3 - 3 = 0
# The coefficients for the polynomial z^4 - 2z^3 + 0z^2 + 0z - 3 are:
coeffs = [1, -2, 0, 0, -3]

# Find all roots of the polynomial using numpy
roots = np.roots(coeffs)

# The variable z represents a physical property, so we are interested in real roots.
# We filter the roots to keep only the ones where the imaginary part is zero.
real_roots = roots[np.isclose(roots.imag, 0)].real

# The problem implies a single answer, and one of the real roots is a clean integer.
# We will select this integer root to find the position x0.
z0 = None
for r in real_roots:
    if np.isclose(r, round(r)):
        z0 = int(round(r))
        break

# The corresponding position x0 is determined by the relation: x0 = 3*z0 - z0^2
# We can now substitute the value of z0 to find x0.
if z0 is not None:
    c1 = 3
    c2 = -1
    x0 = c1 * z0 + c2 * (z0**2)
    
    print("A key step in the solution is solving the polynomial equation z^4 - 2*z^3 - 3 = 0.")
    print(f"An integer root of this equation is z0 = {z0}.")
    print("\nThe position x0 is then found using the relation x0 = 3*z0 - z0^2.")
    print("Substituting the value of z0:")
    print(f"x0 = {c1} * ({z0}) - ({z0})^2")
    print(f"x0 = {c1 * z0} - {z0**2}")
    print(f"x0 = {x0}")
else:
    print("No integer root was found.")
