import math

# The problem simplifies to calculating the value of the expression 1 / (1 + 1/sqrt(2)),
# which is algebraically equivalent to 2 - sqrt(2).

# Define the numbers in the final equation
num1 = 2
num2 = 2

# Calculate the result
result = num1 - math.sqrt(num2)

# Print the final equation as requested
print("The value of (1 - cos(theta_14)) / (1 - cos(theta_34)) is given by the equation:")
print(f"{num1} - sqrt({num2})")

# Print the final numerical value
print("\nThe numerical value of the expression is:")
print(result)