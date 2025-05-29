from sympy import symbols, Eq, solve

# Define symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the problem statement
equations = [
    Eq(F - C, -141),
    C > D,
    Eq(E + K, 305),
    Eq(K, 3 * A),
    Eq(H, 3 * F),
    Eq(I + J, 157),
    Eq(K, 1.5 * C),
    Eq(K - A, 150),
    Eq(F + G, 25),
    Eq(B + C, 153)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter solutions to match the possible values
valid_solution = None
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        valid_solution = sol
        break

# Extract the values in alphabetical order
result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D],
          valid_solution[E], valid_solution[F], valid_solution[G], valid_solution[H],
          valid_solution[I], valid_solution[J], valid_solution[K]]

# Output the result
print(f"<<<{result}>>>")