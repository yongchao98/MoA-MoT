import sympy

# Define symbols
x, y = sympy.symbols('x y')

# The two factors of the solution
factor1 = y**2 - 9
factor2 = x**2 + y**2 - 9

# The general solution is the product of these factors set to zero
general_solution = factor1 * factor2

# We want to output the equation in a readable format
# The factored form is (y**2 - 9)*(x**2 + y**2 - 9) = 0
# The expanded form is y**4 + x**2*y**2 - 9*x**2 - 18*y**2 + 81 = 0

# Print the final equation in factored form with each number.
print(f"(y**2 - {9}) * (x**2 + y**2 - {9}) = {0}")