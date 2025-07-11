import numpy as np
from scipy import integrate

# The problem asks to find the integral of f(x,y,z) = z^2 * (x^2 + y^2)
# over the volume of a cone with base radius 3 and height 2.
# The integral can be solved analytically, yielding the result (108 * pi) / 35.
# As requested, we will output the numbers from this final equation.

numerator_coefficient = 108
denominator = 35

print("--- Analytical Solution Details ---")
print(f"The integral's exact solution has the form (A * pi) / B.")
print(f"The coefficient A in the numerator is: {numerator_coefficient}")
print(f"The denominator B is: {denominator}")
print(f"The symbolic result is: {numerator_coefficient}*pi/{denominator}")
print("-" * 33)

# Now, we will compute the integral numerically using Python.
# For scipy.integrate.tplquad, the function must take arguments in the order (z, y, x).
def f_to_integrate(z, y, x):
    """The function f(x,y,z) = z^2 * (x^2 + y^2)"""
    return z**2 * (x**2 + y**2)

# The integration limits are defined by the cone's geometry.
# z ranges from 0 to 2 (the height of the cone).
z_min = 0
z_max = 2

# For a given z, the radius of the cone's circular cross-section is r(z) = 3 - 1.5*z.
# The y-limits are a function of z.
def y_upper_limit(z):
    return 3 - 1.5 * z

def y_lower_limit(z):
    return -(3 - 1.5 * z)

# The x-limits are a function of z and y, constrained by the circle: x^2 + y^2 <= r(z)^2
def x_upper_limit(z, y):
    radius_sq = (3 - 1.5 * z)**2
    return np.sqrt(radius_sq - y**2)

def x_lower_limit(z, y):
    radius_sq = (3 - 1.5 * z)**2
    return -np.sqrt(radius_sq - y**2)

print("\n--- Numerical Calculation ---")
print("Calculating the integral value using scipy.integrate.tplquad...")

# Perform the triple integration. The result is a tuple (value, estimated_error).
result, error = integrate.tplquad(
    f_to_integrate,
    z_min, z_max,
    y_lower_limit, y_upper_limit,
    x_lower_limit, x_upper_limit
)

print(f"Numerical result from SciPy: {result}")
print(f"Estimated error: {error}")

# For comparison, let's calculate the numerical value from the analytical solution.
analytical_value = numerator_coefficient * np.pi / denominator
print(f"Value from analytical solution (108*pi/35): {analytical_value}")
