from sympy import symbols, Eq, solve

# Define the symbols for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given constraints
equations = [
    Eq(J + K, 271),
    Eq(E + J, 141),
    Eq(A + C, 300),
    Eq(G + K, 230),
    Eq(B - K, -141),
    Eq(D - C, -209),
    Eq(A + I, 78),
    Eq(F - D, 23),
    Eq(E - H, -28),
    Eq(H + K, 198)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Filter the solution to match the possible values
final_solution = {k: v for k, v in solution.items() if v in possible_values}

# Print the final solution in alphabetical order
print([final_solution[A], final_solution[B], final_solution[C], final_solution[D], final_solution[E], final_solution[F], final_solution[G], final_solution[H], final_solution[I], final_solution[J], final_solution[K]])