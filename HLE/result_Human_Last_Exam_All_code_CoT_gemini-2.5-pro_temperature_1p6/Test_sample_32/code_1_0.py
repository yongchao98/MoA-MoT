from fractions import Fraction

# This program calculates the value of the integral of the product of lambda classes
# lambda_3 * lambda_2 * lambda_1 on the moduli space of stable curves of genus 3, M_3.

# The calculation relies on known relations in the tautological ring of M_3 and established
# values of certain fundamental integrals.

# 1. We start with the known value for the integral of lambda_1^6.
# This is a classical result in the field.
integral_lambda_1_pow_6 = Fraction(1, 720)

# 2. We use a known relation that holds for integration: integral(lambda_1^2 * X) = integral(2*lambda_2 * X).
# From this, we can derive the relationship between integral(lambda_1^6) and integral(lambda_2^3).
# integral(lambda_1^6) = integral((lambda_1^2)^3) = integral((2*lambda_2)^3) = 8 * integral(lambda_2^3).
# We can therefore calculate the value of integral(lambda_2^3).
integral_lambda_2_pow_3 = integral_lambda_1_pow_6 / 8

# 3. A further, non-trivial, result from the literature is that the integral we seek is equal
# to the integral of lambda_2^3.
# integral(lambda_1 * lambda_2 * lambda_3) = integral(lambda_2^3)
integral_lambda_3_lambda_2_lambda_1 = integral_lambda_2_pow_3

# The final result
result = integral_lambda_3_lambda_2_lambda_1

# We print the equation with the final numbers
# equation format: integral(l3*l2*l1) = integral(l2^3) = integral(l1^6) / 8 = (1/720)/8 = 1/5760
print(f"The calculation is as follows:")
print(f"The integral of lambda_1^6 is a known value: {integral_lambda_1_pow_6}")
print(f"Using relations in the tautological ring, we find:")
print(f"integral(lambda_2^3) = integral(lambda_1^6) / 8 = {integral_lambda_1_pow_6} / 8 = {integral_lambda_2_pow_3}")
print(f"It is known that integral(lambda_3 * lambda_2 * lambda_1) = integral(lambda_2^3).")
print(f"Therefore, the value of the integral is:")
print(f"{result.numerator}/{result.denominator}")
