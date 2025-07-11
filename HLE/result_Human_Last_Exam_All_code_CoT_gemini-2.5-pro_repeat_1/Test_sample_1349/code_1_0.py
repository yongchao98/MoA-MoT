import math

# The problem asks for the supremum of X, which was found to be
# X_sup = 40 / (25 + 36 * pi^2)

# We define the numbers in this final equation
numerator = 40.0
constant_term_in_denominator = 25.0
pi_squared_coefficient_in_denominator = 36.0

# Calculate pi squared
pi_squared = math.pi**2

# Calculate the denominator
denominator = constant_term_in_denominator + pi_squared_coefficient_in_denominator * pi_squared

# Calculate the final value for the supremum of X
supremum_X = numerator / denominator

# As requested, we print each number in the final equation and then the result.
print("The final expression for the supremum of X is of the form: num / (c1 + c2 * pi^2)")
print(f"num = {numerator}")
print(f"c1 = {constant_term_in_denominator}")
print(f"c2 = {pi_squared_coefficient_in_denominator}")
print(f"pi^2 = {pi_squared}")
print("\nCalculating the final value:")
print(f"Supremum of X = {numerator} / ({constant_term_in_denominator} + {pi_squared_coefficient_in_denominator} * {pi_squared})")
print(f"Supremum of X = {supremum_X}")