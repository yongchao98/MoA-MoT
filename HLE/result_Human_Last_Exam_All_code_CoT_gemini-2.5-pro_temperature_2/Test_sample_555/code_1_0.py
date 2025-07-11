from fractions import Fraction

# We derived the following equation for the minimal expected number of rolls, E:
# E * (2400/2401) = 752/343
# To solve for E, we rearrange the equation:
# E = (752/343) * (2401/2400)

# We use the Fraction class for exact fractional arithmetic.
a = Fraction(752, 343)
b = Fraction(2401, 2400)
E = a * b

# The problem asks for the simplified fraction.
# The Fraction class automatically simplifies the result.
numerator = E.numerator
denominator = E.denominator

# The final equation is E = numerator / denominator
print(f"The minimal expected value is the fraction: {numerator} / {denominator}")
print("The final equation can be written as:")
print(f"{denominator} * E = {numerator}")
