from sympy import symbols, Eq, solve, N

# Define the variable and the equation
p = symbols('p')
equation = Eq(71*p**3 - 2*p - 29, 0)

# Attempt to find rational roots
possible_rational_roots = [1, -1, 29, -29, 1/71, -1/71, 29/71, -29/71]
rational_roots = [root for root in possible_rational_roots if equation.subs(p, root) == 0]

# If no rational roots are found, find numerical solutions
if not rational_roots:
    numerical_roots = solve(equation, p)
    real_roots = [N(root, 4) for root in numerical_roots if root.is_real]
else:
    real_roots = rational_roots

# Print the real roots rounded to 4 decimal places
print(", ".join(f"{root:.4f}" for root in real_roots))