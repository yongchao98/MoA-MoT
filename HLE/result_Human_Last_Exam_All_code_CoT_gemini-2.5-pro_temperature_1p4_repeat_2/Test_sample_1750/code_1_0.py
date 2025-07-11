import math
from fractions import Fraction

# We need to compute C(50, 25) / 2^50.
# C(n, k) is the binomial coefficient "n choose k".

n = 50
k = 25

# Calculate the binomial coefficient C(50, 25)
# math.comb(n, k) calculates n! / (k! * (n-k)!)
numerator = math.comb(n, k)

# Calculate the denominator 2^50
denominator = 2**n

# The problem asks for the result as a fraction.
# We can find the greatest common divisor (GCD) to simplify the fraction.
common_divisor = math.gcd(numerator, denominator)

# Simplify the fraction
simplified_numerator = numerator // common_divisor
simplified_denominator = denominator // common_divisor

# Print the final equation
print(f"I/pi = C({n}, {k}) / 2^{n} = {numerator} / {denominator} = {simplified_numerator} / {simplified_denominator}")

# To fulfill the final answer format requirement.
# Note that this is the value of the integral divided by pi.
final_answer_string = f"{simplified_numerator}/{simplified_denominator}"
# The problem asks for the value as the final answer in a special format.
# Let's construct the output for that format as well, though the print statement above is more informative.
# Example: <<<answer>>>
final_answer = f"<<<{simplified_numerator}/{simplified_denominator}>>>"
# The print statement is the primary output. The final_answer variable is for context.
# print(final_answer)
