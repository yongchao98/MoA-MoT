from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the inequalities
equations = [
    Eq(I - E, -4),
    Eq(F, 2.5 * K),
    Eq(I, 1.5 * B),
    Eq(L - A, 33),
    Eq(E, 2.8 * H),
    Eq(C + M, 111),
    Eq(L - C, -60),
    Eq(H, 2.0 * F),
    Eq(F + L, 41),
    Eq(K - E, -26),
    Eq(A + F, 8)
]

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter the solution to match the possible values
valid_solution = None
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        valid_solution = sol
        break

# Extract the values in alphabetical order
result = [valid_solution[var] for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
print(result)