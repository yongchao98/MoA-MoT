import scipy.integrate

# Define the value of z
z = 0.2

# Define the common part of the integrand
# This is derived from f_Y(y) / (y*(1-y))
# f_Y(y) = 12 * y * (2 - y) * (1 - y)**3
# So f_Y(y) / (y * (1 - y)) = 12 * (2 - y) * (1 - y)**2
integrand_part = lambda y: 12 * (2 - y) * (1 - y)**2

# Define the functions for each of the three integrals
# The integrands are the common part multiplied by L(z, y)
# Case 1: y is in (0, z), L(z,y) = y
func1 = lambda y: integrand_part(y) * y
# Case 2: y is in (z, 1-z), L(z,y) = z
func2 = lambda y: integrand_part(y) * z
# Case 3: y is in (1-z, 1), L(z,y) = 1-y
func3 = lambda y: integrand_part(y) * (1 - y)

# Calculate the first integral from 0 to z
integral_1, _ = scipy.integrate.quad(func1, 0, z)
print(f"The first part of the integral (y from 0 to {z}) is: {integral_1}")

# Calculate the second integral from z to 1-z
integral_2, _ = scipy.integrate.quad(func2, z, 1 - z)
print(f"The second part of the integral (y from {z} to {1-z}) is: {integral_2}")

# Calculate the third integral from 1-z to 1
integral_3, _ = scipy.integrate.quad(func3, 1 - z, 1)
print(f"The third part of the integral (y from {1-z} to 1) is: {integral_3}")

# The final value is the sum of the three integrals
f_z_value = integral_1 + integral_2 + integral_3

print(f"\nThe final calculated value of f_Z({z}) is the sum of these parts:")
print(f"{integral_1} + {integral_2} + {integral_3} = {f_z_value}")
print(f"\nThus, the exact value of f(0.2) is {f_z_value:.3f}")
