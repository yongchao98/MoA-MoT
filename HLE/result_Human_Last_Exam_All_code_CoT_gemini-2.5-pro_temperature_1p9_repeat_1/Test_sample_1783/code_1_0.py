import sympy as sp

# Define the symbols for cylindrical coordinates
r, theta, z = sp.symbols('r theta z')

# Define the function f(x,y,z) = z^2*(x^2+y^2) in cylindrical coordinates
# x^2 + y^2 = r^2, so f(r, theta, z) = z^2 * r^2
f_cylindrical = z**2 * r**2

# The volume element in cylindrical coordinates is r*dz*dr*d(theta).
# The integrand is the function multiplied by the volume element Jacobian, r.
integrand = f_cylindrical * r

# Define the limits of integration for the cone
# Height H=2, Radius R=3.
# The upper z-limit is defined by the cone's surface: z = H*(1 - r/R)
z_limit_upper = 2 * (1 - r/3)
r_limit_upper = 3
theta_limit_upper = 2 * sp.pi

print("Setting up the integral in cylindrical coordinates:")
print(f"Integrand: {integrand}")
print(f"Integration order: dz dr dtheta")
print(f"Limits for z: 0 to {z_limit_upper}")
print(f"Limits for r: 0 to {r_limit_upper}")
print(f"Limits for theta: 0 to {theta_limit_upper}")
print("\nCalculating...")

# Perform the nested integration
# 1. Integrate with respect to z
integral_z = sp.integrate(integrand, (z, 0, z_limit_upper))
# 2. Integrate the result with respect to r
integral_r = sp.integrate(integral_z, (r, 0, r_limit_upper))
# 3. Integrate the result with respect to theta
total_integral = sp.integrate(integral_r, (theta, 0, theta_limit_upper))

# The symbolic result contains the numbers from the final fraction (108, pi, 35)
symbolic_result = total_integral
# The numerical value is the floating point approximation
numerical_result = total_integral.evalf()

print("The integral has been calculated.")
# The following line outputs the numbers 108 and 35 which are part of the final equation/answer.
print(f"The exact symbolic result is: {symbolic_result}")
print(f"The approximate numerical value is: {numerical_result}")
