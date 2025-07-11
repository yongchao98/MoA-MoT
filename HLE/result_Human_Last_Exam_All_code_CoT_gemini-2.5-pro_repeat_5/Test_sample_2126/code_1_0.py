import math
import cmath

# Based on the analysis, the quantity can be calculated as (-6i)^(1/2) * G'(0).

# 1. Define the parameters of the final equation
c = 6
# From the derivation, G'(0) = -2 / sqrt(6)
G_prime_0_numerator = -2
G_prime_0_denominator = math.sqrt(6)
G_prime_0 = G_prime_0_numerator / G_prime_0_denominator

# 2. Calculate the complex factor (-ic)^(1/2)
complex_factor = cmath.sqrt(-1j * c)

# 3. Calculate the final result
result = complex_factor * G_prime_0

# 4. Output the numbers used in the final equation and the result
print(f"The calculation is based on the formula: Q = (-i*c)^(1/2) * G'(0)")
print(f"Wave speed c = {c}")
print(f"G'(0) = {G_prime_0_numerator} / sqrt({G_prime_0_denominator**2}) = {G_prime_0}")
print(f"Complex factor (-{c}i)^(1/2) = {complex_factor}")
print(f"Final Quantity Q = ({complex_factor}) * ({G_prime_0})")
print(f"Q = {result}")

# The final answer in the required format
final_answer_str = f"<<<{result}>>>"
print(final_answer_str)