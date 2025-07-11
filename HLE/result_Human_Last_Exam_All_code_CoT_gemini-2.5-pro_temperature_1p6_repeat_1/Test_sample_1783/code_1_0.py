import math

# These are the coefficients derived from the analytical solution of the integral.
numerator_coefficient = 108
denominator = 35

# The symbol 'pi' is also a number in our final equation.
pi = math.pi

# Calculate the final result of the integral.
result = (numerator_coefficient * pi) / denominator

# As requested, here is the final equation showing each number used.
print(f"The equation for the integral's value is ({numerator_coefficient} * \u03C0) / {denominator}")
print(f"The final result is: {result}")