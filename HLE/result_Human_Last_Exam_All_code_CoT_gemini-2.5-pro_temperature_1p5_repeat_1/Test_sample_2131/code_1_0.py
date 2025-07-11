# This script calculates the membrane's deflection y(0) based on the derived formula.

# Define the numerical components of the formula y(0) = (1/3) * (9/4)^(8/5)
coefficient_numerator = 1
coefficient_denominator = 3
base_numerator = 9
base_denominator = 4
exponent_numerator = 8
exponent_denominator = 5

# Calculate the values
coefficient = coefficient_numerator / coefficient_denominator
base = base_numerator / base_denominator
exponent = exponent_numerator / exponent_denominator

# Calculate the final deflection y(0)
deflection_y0 = coefficient * (base ** exponent)

# Print the equation with its components and the final calculated value
print("The membrane's deflection at x=0, y(0), is calculated using the following equation:")
print(f"y(0) = ({coefficient_numerator}/{coefficient_denominator}) * ({base_numerator}/{base_denominator})^({exponent_numerator}/{exponent_denominator})")
print("\nThe calculated deflection is:")
print(deflection_y0)