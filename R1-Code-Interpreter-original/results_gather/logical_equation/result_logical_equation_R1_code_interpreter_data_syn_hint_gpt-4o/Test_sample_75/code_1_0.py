from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the problem statement
equations = [
    Eq(L - I, -23),
    Eq(I + M, 44),
    Eq(D + K, 98),
    Eq(F + I, 73),
    Eq(F - J, 38),
    Eq(H - F, 5),
    Eq(J, 1.4 * L),
    Eq(B - C, 26),
    Eq(B, 2.4 * A),
    Eq(C + F, 55)
]

# Define the list of possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter solutions that match the possible values
valid_solution = None
for sol in solution:
    if all(val in possible_values for val in sol.values()):
        valid_solution = sol
        break

# Extract the solution in alphabetical order
if valid_solution:
    result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D], valid_solution[E],
              valid_solution[F], valid_solution[G], valid_solution[H], valid_solution[I], valid_solution[J],
              valid_solution[K], valid_solution[L], valid_solution[M]]
    print(result)