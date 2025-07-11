from fractions import Fraction

# The problem reduces to solving the linear equation for E1 (the expected number of rolls):
# E1 = a + b * E1
# where a = 752/343 and b = 9/21609

# Coefficients from the equation: E1 = (752/343) + (9/21609) * E1
a = Fraction(752, 343)
b = Fraction(9, 21609)

# We rearrange the equation to solve for E1:
# E1 - b * E1 = a
# E1 * (1 - b) = a
# E1 = a / (1 - b)

one_minus_b = 1 - b
E1 = a / one_minus_b

print(f"The final equation to be solved is E1 = {a.numerator}/{a.denominator} + ({b.numerator}/{b.denominator}) * E1")
print(f"Solving for E1 gives E1 = ({a.numerator}/{a.denominator}) / (1 - {b.numerator}/{b.denominator})")
print(f"E1 = ({a.numerator}/{a.denominator}) / ({one_minus_b.numerator}/{one_minus_b.denominator})")
print(f"The minimal expected number of rolls is: {E1.numerator}/{E1.denominator}")