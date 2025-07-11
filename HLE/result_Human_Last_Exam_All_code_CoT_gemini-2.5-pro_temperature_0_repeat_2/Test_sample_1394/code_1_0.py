import sympy

# Define the symbols
y, x, C = sympy.symbols('y x C')

# The general solution
general_solution = sympy.Eq(y**2, C*x**2 + C**2 + C + 9)

# Print the solution in a readable format
print("The general solution of the differential equation is:")
print(sympy.pretty(general_solution, use_unicode=False))
print("\nwhere C is an arbitrary constant.")

# As requested, output each number in the final equation
print("\nThe numbers in the final equation y**2 = C*x**2 + C**2 + C + 9 are:")
# The power of y
print(2)
# The power of x
print(2)
# The power of C
print(2)
# The constant term
print(9)
