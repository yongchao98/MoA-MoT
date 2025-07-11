import numpy as np
from scipy import integrate

# The function to integrate is f(x,y,z) = z^2 * (x^2 + y^2).
# In cylindrical coordinates, this is f(r, theta, z) = z^2 * r^2.
# The volume element dV is r*dr*d_theta*dz.
# So the integrand is (z^2 * r^2) * r = r^3 * z^2.

# Note on scipy.integrate.tplquad order:
# It integrates func(z, y, x) dz dy dx.
# We map our variables as: x -> theta, y -> z, z -> r
# So the function provided to tplquad must be func(r, z, theta)
integrand = lambda r, z, theta: r**3 * z**2

# Limits for theta (our 'x')
theta_min = 0
theta_max = 2 * np.pi

# Limits for z (our 'y'), which can be functions of theta
z_min = lambda theta: 0
z_max = lambda theta: 2

# Limits for r (our 'z'), which can be functions of theta and z
r_min = lambda theta, z: 0
r_max = lambda theta, z: 3 * (1 - z/2.0)

# Calculate the integral
result, error = integrate.tplquad(integrand,
                                   theta_min, theta_max,
                                   z_min, z_max,
                                   r_min, r_max)

# --- Output ---
print("The problem is to find the integral of f(x,y,z) = z^2*(x^2+y^2) inside a cone of radius 3 and height 2.")
print("This is solved by setting up a triple integral in cylindrical coordinates.\n")
print("The integral is set up as follows:")
print("Integral(theta=0 to 2*pi) { Integral(z=0 to 2) { Integral(r=0 to 3*(1-z/2)) of (r**3 * z**2) dr } dz } d_theta\n")

# Print the final equation with each number as requested
symbolic_result = "108*pi/35"
print(f"The final equation and its symbolic result is:")
print(f"Integral from theta=0 to 2*pi of "
      f"Integral from z=0 to 2 of "
      f"Integral from r=0 to 3*(1-z/2) of (r**3 * z**2) dr dz d_theta = {symbolic_result}")
print("\nCalculating the numerical value:")
print(f"The symbolic result is 108*pi/35 ≈ {108*np.pi/35:.4f}")
print(f"The numerical result from scipy.integrate.tplquad is: {result:.4f}")
