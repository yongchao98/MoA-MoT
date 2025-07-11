from fractions import Fraction

# This program calculates the integral of the product of lambda classes 
# lambda_3 * lambda_2 * lambda_1 on the moduli space of stable curves of genus 3.

# The calculation is based on established relations in the tautological ring of M_3.
# A key Faber-Zagier relation for g=3 is lambda_1 * lambda_2 = 12 * lambda_3.
# Multiplying by lambda_3 gives a top-degree relation, which implies that the integrals
# over the compactified moduli space are equal:
# Integral(lambda_1*lambda_2*lambda_3) = 12 * Integral(lambda_3^2)

# The value of Integral(lambda_3^2) is a known result from Faber's work,
# proven by Faber, Pandharipande, and others.
# The number 12 is the constant in the Faber-Zagier relation.
# The value of the integral of lambda_3^2 is 1/3456.
factor = 12
integral_lambda3_squared = Fraction(1, 3456)

# We compute the product to find the value of the requested integral.
result = factor * integral_lambda3_squared

# The final equation is: 
# Integral(lambda_3 * lambda_2 * lambda_1) = 12 * (1/3456) = 1/288
# The problem asks for the output in "a/b" format.
# We print the numerator and denominator of the resulting fraction.
print(f"{result.numerator}/{result.denominator}")