from sympy import symbols, Eq, solve, S

# Define the variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I', integer=True)

# Define the possible values
possible_values = {1, 3, 4, 9, 16, 27, 36, 80, 121}

# Define the equations and inequalities
equations = [
    Eq(A - D, -77),
    Eq(I, 3 * E),
    Eq(E - H, -27),
    Eq(G - A, -2),
    Eq(H, 4 * E)
]

inequalities = [
    H > A,
    I > B
]

# Solve the equations and inequalities
solutions = solve(equations + inequalities, (A, B, C, D, E, F, G, H, I), domain=S.Integers)

# Filter solutions to match possible values
valid_solution = None
for sol in solutions:
    if all(value in possible_values for value in sol):
        valid_solution = sol
        break

# Output the solution in the required format
if valid_solution:
    result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D],
              valid_solution[E], valid_solution[F], valid_solution[G], valid_solution[H], valid_solution[I]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")