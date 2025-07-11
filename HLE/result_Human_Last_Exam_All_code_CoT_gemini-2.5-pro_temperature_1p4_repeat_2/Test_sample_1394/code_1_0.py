import sympy as sp

# Define the variables
x, y, C = sp.symbols('x y C')

# The two families of solutions derived from the factored differential equation
solution1 = 2*y + 6*sp.log(sp.Abs(y - 3)) - x**2 - C
solution2 = 2*y - 6*sp.log(sp.Abs(y + 3)) + x**2 - C

# The general solution is the product of the two solutions set to zero
general_solution = sp.Eq(solution1 * solution2, 0)

# Print the general solution in a readable format
# The general solution is (2*y + 6*ln(|y-3|) - x**2 - C) * (2*y - 6*ln(|y+3|) + x**2 - C) = 0
print(f"The general solution is given by the equation:")
print(f"({2}*y + {6}*log(|y - {3}|) - x**{2} - C) * ({2}*y - {6}*log(|y + {3}|) + x**{2} - C) = {0}")
