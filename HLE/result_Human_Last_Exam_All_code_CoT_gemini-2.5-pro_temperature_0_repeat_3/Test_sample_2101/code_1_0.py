import math

# The difference between the probability that an electron escapes through the hypotenuse
# and the probability that it escapes through either of the two legs can be
# simplified to the expression: 2*sqrt(2) - 3.
# Here, we calculate its numerical value.

# Define the constant value for the square root of 2.
val_sqrt2 = math.sqrt(2)

# Define the terms in the final simplified equation.
num1 = 2
num2 = 3

# Calculate the final result.
result = num1 * val_sqrt2 - num2

print("The simplified expression for the difference in probabilities is 2*sqrt(2) - 3.")
print("The final equation with the numerical values is:")
print(f"{num1} * {val_sqrt2} - {num2} = {result}")