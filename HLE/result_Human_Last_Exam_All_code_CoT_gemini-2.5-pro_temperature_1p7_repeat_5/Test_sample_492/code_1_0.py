import sympy

# Plan: Solve the equation for the emergence of the giant component.
# The derived equation is (2 * c^2) / 5 = 1.

# Define the symbol 'c' for the critical time.
c = sympy.Symbol('c', positive=True)

# Define the coefficients of the final equation.
# Equation is of the form: (a * c**2) / b = d
a = 2
b = 5
d = 1

# Construct the equation using sympy.
final_equation = sympy.Eq(a * c**2 / b, d)

# As per the instructions, we output the numbers in the final equation.
print("Based on the derivation, the critical time 'c' is the solution to the following equation:")
print(f"({a} * c**2) / {b} = {d}")
print("-" * 35)

# Solve the equation for 'c'. Since we defined 'c' as positive,
# sympy will return only the positive root.
solution = sympy.solve(final_equation, c)

# The result is a list, so we extract the first element.
c_value = solution[0]

# Print the final exact answer.
print(f"The exact value for the emergence time c is: {c_value}")

# Print the numerical approximation of the answer.
print(f"The approximate numerical value is: {c_value.evalf()}")