import math

# The integral simplifies to (pi / 2**49) * C(50, 25)
# We will compute the numerator and denominator of the fraction.

# Numerator is C(50, 25)
numerator = math.comb(50, 25)

# Denominator is 2**49
denominator = 2**49

# Simplify the fraction by dividing by their greatest common divisor
common_divisor = math.gcd(numerator, denominator)

final_numerator = numerator // common_divisor
final_denominator = denominator // common_divisor

# The problem asks to write the answer strictly as a fraction.
# The final answer is (final_numerator / final_denominator) * pi.
# We will print the components of this fraction.
print(f"The integral evaluates to ({final_numerator}/{final_denominator}) * pi")
print(f"Final Answer Numerator: {final_numerator}")
print(f"Final Answer Denominator: {final_denominator}")
