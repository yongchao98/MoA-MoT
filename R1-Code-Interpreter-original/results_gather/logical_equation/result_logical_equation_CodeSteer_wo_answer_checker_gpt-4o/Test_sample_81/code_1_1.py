from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(F + G, 46),
    Eq(J - M, 80),
    Eq(J - F, 86),
    Eq(H - L, 26),
    Eq(B + M, 31),
    Eq(B, 1.5 * F),
    Eq(C, 4.0 * A),
    Eq(L, 1.5 * M),
    Eq(G, 3.6 * F),
    Eq(D, 1.5 * E),
    Eq(L - I, -21)
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter solutions to match possible values
valid_solution = None
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        valid_solution = sol
        break

# Extract the values in alphabetical order
result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D],
          valid_solution[E], valid_solution[F], valid_solution[G], valid_solution[H],
          valid_solution[I], valid_solution[J], valid_solution[K], valid_solution[L],
          valid_solution[M]]

# Output the result
print(f"<<<{result}>>>")