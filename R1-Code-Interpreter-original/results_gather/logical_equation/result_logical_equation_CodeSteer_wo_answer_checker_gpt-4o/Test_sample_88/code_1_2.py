from sympy import symbols, Eq, solve, S

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Define the constraints
constraints = [
    Eq(J, 3.5 * D),
    Eq(K - F, -12),
    Eq(G + M, 26),
    Eq(J + K, 31),
    Eq(C - D, 13),
    Eq(G, 1.6 * M),
    L > G,
    Eq(A + L, 124),
    Eq(M - A, -18),
    Eq(C, 3.0 * H),
    J > D
]

# Solve the system of equations and inequalities
solution = solve(constraints, (A, B, C, D, E, F, G, H, I, J, K, L, M), domain=S.Integers)

# Filter solutions to match possible values
valid_solution = None
for sol in solution:
    if all(val in possible_values for val in sol):
        valid_solution = sol
        break

# Print the valid solution
if valid_solution:
    print(f"<<<{list(valid_solution)}>>>")
else:
    print("No valid solution found.")