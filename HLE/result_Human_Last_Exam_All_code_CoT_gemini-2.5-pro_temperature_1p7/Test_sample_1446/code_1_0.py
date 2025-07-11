# Plan:
# 1. The problem is interpreted as finding the critical exponent ν for the φ⁴ theory.
#    The "G₄" is taken to refer to the quartic term in this theory.
# 2. A key concept in this theory is the upper critical dimension, d=4. At this
#    dimension, the critical exponents take on their simple, rational mean-field values.
# 3. The mean-field value for the critical exponent ν (nu) is precisely 1/2.
# 4. This script calculates this value and prints the equation as requested.

# Numerator and denominator for the mean-field value of nu
numerator = 1
denominator = 2

# Calculation of the critical exponent nu
nu_value = numerator / denominator

# The final equation requires printing each number involved in the calculation.
print(f"The calculation for the critical exponent ν is based on mean-field theory.")
print(f"The equation is: {numerator} / {denominator} = {nu_value}")
