import math

# The problem reduces to finding the diagonal length of the largest regular pentagon
# that can be inscribed in a unit square.
# The value r is given by the formula r = sqrt(2 - (2*sqrt(5))/5).

# Define the constants in the equation
a = 2
b = 5

# Calculate the components of the formula
sqrt_b = math.sqrt(b)
term = a * sqrt_b / b
expression_under_sqrt = a - term
r = math.sqrt(expression_under_sqrt)

# Output the equation with the numbers and the final result
print(f"The final equation for r is: r = sqrt({a} - ({a} * sqrt({b})) / {b})")
print(f"The value of r is: {r}")
