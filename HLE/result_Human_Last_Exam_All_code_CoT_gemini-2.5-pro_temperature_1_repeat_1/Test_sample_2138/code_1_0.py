import sympy

# The analytical value of the integral is (pi^2 * sqrt(3)) / 27.
# We define the numbers in this expression to construct it.
power_of_pi = 2
base_of_sqrt = 3
denominator = 27

# Use sympy to define the symbolic expression for the analytical value.
# sympy.pi is the symbolic constant for pi.
# sympy.sqrt() is the symbolic square root function.
analytical_expression = (sympy.pi**power_of_pi * sympy.sqrt(base_of_sqrt)) / denominator

# Print the final equation, showing how the numbers form the expression
# and the resulting symbolic representation.
print(f"(pi**{power_of_pi} * sqrt({base_of_sqrt})) / {denominator} = {analytical_expression}")

# We can also print its numerical approximation.
print(f"The numerical value is approximately: {analytical_expression.evalf()}")