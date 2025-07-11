import math

# The value of y(0) is given by the analytical formula derived from solving the Lagrange equation.
# The formula is: y(0) = (3/4) * (9/4)^(3/5)

# Define the constants from the formula
coeff_a = 3
coeff_b = 4
base_num = 9
base_den = 4
exp_num = 3
exp_den = 5

# Calculate the base and exponent
base = base_num / base_den
exponent = exp_num / exp_den

# Calculate the final result
result = (coeff_a / coeff_b) * (base ** exponent)

# Print the explanation, the formula with its components, and the final numerical result.
print("The deflection of the membrane at x = 0 is given by the formula:")
print(f"y(0) = ({coeff_a}/{coeff_b}) * ({base_num}/{base_den})**({exp_num}/{exp_den})")
print("\nThis simplifies to:")
print(f"y(0) = (3/4) * (2.25)**(0.6)")
print(f"\nThe numerical value for the deflection y(0) is:")
print(result)