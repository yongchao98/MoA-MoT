import sympy

# Define the symbols for our coordinates
r, theta, z = sympy.symbols('r theta z')

# Define the function in cylindrical coordinates, multiplied by the Jacobian 'r'
integrand = z**2 * r**2 * r

print(f"The integrand in cylindrical coordinates (f*r) is: {integrand}")

# Define the limits of integration
r_limit_upper = (sympy.S(3)/2) * (2 - z)
z_limit_lower = 0
z_limit_upper = 2
theta_limit_lower = 0
theta_limit_upper = 2 * sympy.pi

print(f"\nIntegration limits:")
print(f"theta: from {theta_limit_lower} to {theta_limit_upper}")
print(f"z: from {z_limit_lower} to {z_limit_upper}")
print(f"r: from 0 to {r_limit_upper}")


# 1. Integrate with respect to r
integral_after_r = sympy.integrate(integrand, (r, 0, r_limit_upper))
print(f"\nAfter integrating with respect to r, we get:\n{integral_after_r}")

# 2. Integrate the result with respect to z
integral_after_z = sympy.integrate(integral_after_r, (z, z_limit_lower, z_limit_upper))
print(f"\nAfter integrating that result with respect to z, we get:\n{integral_after_z}")

# 3. Integrate the final result with respect to theta
final_integral = sympy.integrate(integral_after_z, (theta, theta_limit_lower, theta_limit_upper))
print(f"\nFinally, after integrating with respect to theta, the value of the integral is:\n{final_integral}")

# The user wants to see the numbers in the final equation.
# The result is of the form C * pi. We can extract the coefficient C.
# The final result is (54/35) * 2*pi = 108*pi/35
# We can print this composition to satisfy the prompt.
coeff, _ = final_integral.as_coeff_Mul(sympy.pi)
numerator, denominator = sympy.fraction(coeff)
print(f"\nThe final equation can be written as: ({numerator} * pi) / {denominator}")

final_numeric_value = final_integral.evalf()
print(f"\nThe numerical value is approximately: {final_numeric_value}")
print(f"\n<<< {final_integral} >>>")