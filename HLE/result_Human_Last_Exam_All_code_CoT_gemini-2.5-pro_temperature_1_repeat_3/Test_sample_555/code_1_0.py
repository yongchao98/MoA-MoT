from fractions import Fraction

# We have a system of linear equations for the expected number of rolls E_M.
# E_1 = 1 + E_7
# E_7 = 1 + (9/49) * E_9
# E_9 = 1 + (1/21) * E_3
# E_3 = 1 + (1/21) * E_1
# Let's solve this system for E_1 by substitution.

# Let E_1 be represented by a variable 'x'. We will find expressions of the
# form 'c + d*x' for each E_M.

# For E_3: E_3 = 1 + (1/21) * E_1
c3 = Fraction(1)
d3 = Fraction(1, 21)

# For E_9: E_9 = 1 + (1/21) * E_3 = 1 + (1/21) * (c3 + d3*E_1)
c9 = 1 + Fraction(1, 21) * c3
d9 = Fraction(1, 21) * d3

# For E_7: E_7 = 1 + (9/49) * E_9 = 1 + (9/49) * (c9 + d9*E_1)
c7 = 1 + Fraction(9, 49) * c9
d7 = Fraction(9, 49) * d9

# For E_1: E_1 = 1 + E_7 = 1 + (c7 + d7*E_1)
# This gives the final equation for E_1:
# E_1 = (1 + c7) + d7*E_1
# E_1 * (1 - d7) = 1 + c7
# E_1 = (1 + c7) / (1 - d7)

final_c = 1 + c7
final_d = d7
final_denominator = 1 - final_d

# Calculate the final result
E1 = final_c / final_denominator

print("The problem reduces to solving a system of linear equations for E_1, the expected number of rolls.")
print("The final equation for E_1 is derived by substitution:")
print(f"E_1 * (1 - {final_d}) = {final_c}")
print(f"E_1 * ({final_denominator}) = {final_c}")
print(f"E_1 = {final_c} / {final_denominator}")
print()
print(f"The minimal expected value of rolls is {E1}.")
print(f"As a simplified fraction, this is: {E1.numerator}/{E1.denominator}")

<<<329/150>>>