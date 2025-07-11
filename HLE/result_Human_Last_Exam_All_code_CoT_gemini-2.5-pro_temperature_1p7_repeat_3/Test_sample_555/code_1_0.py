from fractions import Fraction

# This script solves the system of linear equations for the expected value E.
# Let E be E_1.
# The system is:
# E_3 = 1 + (1/21) * E
# E_9 = 1 + (3/63) * E_3
# E_7 = 1 + (9/49) * E_9
# E   = 1 + E_7
# We solve this by substitution. An intermediate expectation E_i can be expressed as:
# E_i = A_i * E + B_i

# For E_3:
A3 = Fraction(1, 21)
B3 = Fraction(1)

# For E_9, substituting E_3:
# E_9 = 1 + (3/63) * (A3*E + B3)
A9 = Fraction(3, 63) * A3
B9 = Fraction(1) + Fraction(3, 63) * B3

# For E_7, substituting E_9:
# E_7 = 1 + (9/49) * (A9*E + B9)
A7 = Fraction(9, 49) * A9
B7 = Fraction(1) + Fraction(9, 49) * B9

# Finally, for E:
# E = 1 + E_7 = 1 + (A7*E + B7)
# E * (1 - A7) = 1 + B7
# E = (1 + B7) / (1 - A7)

# Calculate the terms
numerator_term = 1 + B7
denominator_term = 1 - A7

# The final calculation is a division of two fractions, which is equivalent to a multiplication
multiplicand1 = numerator_term
multiplicand2 = Fraction(denominator_term.denominator, denominator_term.numerator)

# Perform the final multiplication
final_result_unsimplified = multiplicand1 * multiplicand2
final_result_simplified = final_result_unsimplified

print("The minimal expected value, E, is calculated by solving the recurrence relations.")
print("The final calculation for E is:")
print(f"E = ({multiplicand1.numerator}/{multiplicand1.denominator}) * ({multiplicand2.numerator}/{multiplicand2.denominator})")
print(f"E = {final_result_unsimplified.numerator}/{final_result_unsimplified.denominator}")
print()
print("The simplified result is:")
print(f"E = {final_result_simplified.numerator}/{final_result_simplified.denominator}")