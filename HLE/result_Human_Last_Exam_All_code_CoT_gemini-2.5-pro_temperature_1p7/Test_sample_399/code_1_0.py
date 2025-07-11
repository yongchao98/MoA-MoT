from fractions import Fraction

# We are calculating x(1) using the derived formula:
# x(t) = -8 / (4*t^2 - 8*t + 1)
t = 1

# Define the components of the formula
numerator = -8
# Calculate the denominator by substituting t=1
term1 = 4 * (t**2)
term2 = -8 * t
term3 = 1
denominator = term1 + term2 + term3

# Calculate the final result
result = Fraction(numerator, denominator)

# Print the step-by-step calculation
print(f"We need to find x(1) using the equation x(t) = -8 / (4*t^2 - 8*t + 1).")
print(f"Substituting t = 1:")
print(f"x(1) = {numerator} / (4*({t})^2 - 8*({t}) + 1)")
print(f"x(1) = {numerator} / ({term1} - 8 + {term3})")
print(f"x(1) = {numerator} / ({denominator})")
print(f"The final result is: x(1) = {result}")
