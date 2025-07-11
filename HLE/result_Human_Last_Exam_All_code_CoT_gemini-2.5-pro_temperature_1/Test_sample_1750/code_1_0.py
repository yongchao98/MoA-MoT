import math
from fractions import Fraction

# Calculate the components of the fractional coefficient of pi
# Numerator is the binomial coefficient "50 choose 25"
numerator = math.comb(50, 25)

# Denominator is 2 to the power of 50
denominator = 2**50

# The equation for the integral I is:
# I = (numerator / denominator) * pi
print("The value of the integral is given by the equation:")
print(f"I = ({numerator} / {denominator}) * pi")

# The problem asks for the answer as a fraction. We provide the simplified
# fractional coefficient of pi.
coefficient_fraction = Fraction(numerator, denominator)

final_numerator = coefficient_fraction.numerator
final_denominator = coefficient_fraction.denominator

print("\nThe simplified fractional coefficient of pi is:")
print(f"Numerator: {final_numerator}")
print(f"Denominator: {final_denominator}")
print(f"So, the answer is ({final_numerator}/{final_denominator}) * pi")

# Format the final answer as requested
final_answer_str = f"{final_numerator}/{final_denominator}"
print(f"\nFinal answer in the requested format:")
print(f"<<<{final_answer_str}>>>")
