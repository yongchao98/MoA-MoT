import sympy as sp

# Define the symbol 'a' as a positive real number
a = sp.Symbol('a', positive=True)

# The derived equation for the arc length L in terms of 'a' is L = 3 * a^(2/3).
# We are given that L = 3/2.

# Set up the equation: 3 * a^(2/3) = 3/2
lhs_coeff = 3
power_num = 2
power_den = 3
rhs_num = 3
rhs_den = 2

# We use sp.S to represent rational numbers precisely for accurate symbolic computation
equation = sp.Eq(lhs_coeff * a**(sp.S(power_num)/power_den), sp.S(rhs_num)/rhs_den)

# Print the equation that needs to be solved
print("The equation relating the arc length to the constant 'a' is:")
print(f"{lhs_coeff} * a^({power_num}/{power_den}) = {rhs_num}/{rhs_den}")

# Solve the equation for 'a'
solutions = sp.solve(equation, a)

# The result is a list, and since 'a' is positive, we take the single solution.
solution_a = solutions[0]

# Print the solution for 'a'
print("\nThe value of 'a' is:")
sp.pprint(solution_a)