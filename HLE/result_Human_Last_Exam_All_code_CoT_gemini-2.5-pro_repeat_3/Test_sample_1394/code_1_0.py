import sympy

# Define the variables and the constant
x, y = sympy.symbols('x y')
C = sympy.Symbol('C')

# The general solution to the differential equation is given by the following implicit equation:
# 4*y**2 - 2*C*x**2 - C**2 - 2*C - 36 = 0
# We can represent this equation using sympy
general_solution = 4*y**2 - 2*C*x**2 - C**2 - 2*C - 36

# Print the general solution
print("The general solution of the differential equation is:")
print(f"{general_solution} = 0, where C is an arbitrary constant.")

# As requested, output each number in the final equation.
# The equation is of the form: A*y**2 + B*C*x**2 + D*C**2 + E*C + F = 0
A = 4
B = -2
D = -1
E = -2
F = -36

print("\nThe numbers (coefficients) in the final equation are:")
print(f"A (coefficient of y**2): {A}")
print(f"B (coefficient of C*x**2): {B}")
print(f"D (coefficient of C**2): {D}")
print(f"E (coefficient of C): {E}")
print(f"F (constant term): {F}")