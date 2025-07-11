import sympy

# In this variant of the Erdos-Renyi random graph, we need to find the critical
# time 'c' for the emergence of the giant component. This occurs when the
# average degree of the graph becomes 1.

# Through analysis, we can derive the average degree k at time t as:
# k(t) = t^2 / 3
# We set k(c) = 1 to find the critical time c.

# Define 'c' as a symbolic variable to solve the equation.
c = sympy.Symbol('c')

# Define the equation for the critical time c: c^2 / 3 = 1
# The numbers in the equation are the exponent 2, the denominator 3, and the right-hand side 1.
exponent = 2
denominator = 3
rhs_val = 1
equation = sympy.Eq(c**exponent / denominator, rhs_val)

# Solve the equation for c. We are interested in the positive solution as time is non-negative.
solutions = sympy.solve(equation, c)
positive_solution = [sol for sol in solutions if sol.is_positive][0]

# Print the derivation of the final equation step-by-step.
print("The final equation for the critical time c is derived from the condition that the average degree is 1.")
print(f"The average degree at time c is c^{exponent} / {denominator}.")
print("Setting this to 1 gives the equation:")
print("")

# Print the equation with each number explicitly shown.
print(f"c^{exponent} / {denominator} = {rhs_val}")

# Rearrange the equation to solve for c^2.
new_rhs = rhs_val * denominator
print(f"c^{exponent} = {new_rhs}")

# Take the square root to find c.
print(f"c = sqrt({new_rhs})")
print("")

# Print the final exact value of c.
print(f"The exact value of c is: {positive_solution}")