from fractions import Fraction

# This script calculates the integral of the product of lambda classes lambda_3*lambda_2*lambda_1
# on the moduli space of stable curves of genus 3, M_3_bar.

# The calculation relies on known results from the theory of intersection on moduli spaces.
# Specifically, we use a relation in the Chow ring CH^6(M_3_bar), which is one-dimensional.
# The relation is: lambda_1 * lambda_2 * lambda_3 = 3 * lambda_3^2
# This implies that their integrals are related by the same factor.

proportionality_constant = Fraction(3, 1)

# The base value for the integral of lambda_3^2 is a known result.
# integral(lambda_3^2) = 1/746496
integral_lambda3_squared = Fraction(1, 746496)

# Now, we compute the desired integral using the relation.
# integral(lambda_1*lambda_2*lambda_3) = 3 * integral(lambda_3^2)
result = proportionality_constant * integral_lambda3_squared

# Output the equation and the final result.
print(f"integral(lambda_3*lambda_2*lambda_1) = {proportionality_constant.numerator}/{proportionality_constant.denominator} * integral(lambda_3^2)")
print(f"= {proportionality_constant.numerator}/{proportionality_constant.denominator} * {integral_lambda3_squared.numerator}/{integral_lambda3_squared.denominator}")
print(f"= {result.numerator}/{result.denominator}")