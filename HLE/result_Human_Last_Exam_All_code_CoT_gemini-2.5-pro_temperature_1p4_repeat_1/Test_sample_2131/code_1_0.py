import math

# The final equation for y(0) is y(0) = (3/4) * (9/4)^(3/5)
# We will define the numbers in this equation as variables.
num1 = 3
den1 = 4
num2 = 9
den2 = 4
num_pow = 3
den_pow = 5

# Calculate the value of y(0)
base = num2 / den2
exponent = num_pow / den_pow
result = (num1 / den1) * math.pow(base, exponent)

# Output the equation and the final result.
# As requested, we show each number in the final equation.
print(f"The final equation for the deflection at x=0 is y(0) = ({num1}/{den1}) * ({num2}/{den2})^({num_pow}/{den_pow})")
print("The calculated value for the deflection y(0) is:")
print(result)
