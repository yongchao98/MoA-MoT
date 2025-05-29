from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
equations = [
    Eq(H, 2.4 * D),
    Eq(G - A, -9),
    Eq(E, 4.0 * G),
    Eq(F - G, 8),
    Eq(B + E, 78),
    Eq(M - K, 40),
    Eq(D + M, 55),
    Eq(A, 3.2 * K),
    Eq(I - C, 94),
    Eq(G + L, 43),
    Eq(L - A, 20)
]

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
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