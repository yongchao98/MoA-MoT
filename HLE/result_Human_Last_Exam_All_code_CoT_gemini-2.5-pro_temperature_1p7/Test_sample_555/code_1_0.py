from fractions import Fraction

# The problem is modeled by the following system of linear equations for the expected number of rolls:
# E = 2 + (9/49) * E_9
# E_9 = 1 + (3/63) * E_3
# E_3 = 1 + (1/21) * E
#
# We solve this system for E.
# From the last equation, substitute E_3 into the second:
# E_9 = 1 + (3/63) * (1 + (1/21) * E)
# E_9 = 1 + Fraction(3, 63) * (1 + Fraction(1, 21) * E)
# E_9 = 1 + Fraction(1, 21) + Fraction(1, 441) * E
# E_9 = Fraction(22, 21) + Fraction(1, 441) * E
#
# Now substitute this expression for E_9 into the first equation:
# E = 2 + (9/49) * (22/21 + (1/441) * E)
# E = 2 + Fraction(9, 49) * Fraction(22, 21) + Fraction(9, 49) * Fraction(1, 441) * E
# E = 2 + Fraction(198, 1029) + Fraction(9, 21609) * E
#
# Simplify the fractions:
# Fraction(198, 1029) simplifies to 66/343.
# Fraction(9, 21609) simplifies to 1/2401.
# E = 2 + Fraction(66, 343) + Fraction(1, 2401) * E
#
# Rearrange to solve for E:
# E * (1 - Fraction(1, 2401)) = 2 + Fraction(66, 343)
# E * Fraction(2400, 2401) = Fraction(2*343 + 66, 343)
# E * Fraction(2400, 2401) = Fraction(752, 343)
#
# E = Fraction(752, 343) / Fraction(2400, 2401)

# Let's calculate the final result using the fractions module.
numerator_term = Fraction(752, 343)
denominator_term = Fraction(2400, 2401)
expected_value = numerator_term / denominator_term

print("The minimal expected number of rolls is E, calculated from the equation:")
print("E * (2400 / 2401) = 752 / 343")
print("\nThe solution is the simplified fraction:")
print(f"Numerator: {expected_value.numerator}")
print(f"Denominator: {expected_value.denominator}")
print(f"E = {expected_value.numerator} / {expected_value.denominator}")
