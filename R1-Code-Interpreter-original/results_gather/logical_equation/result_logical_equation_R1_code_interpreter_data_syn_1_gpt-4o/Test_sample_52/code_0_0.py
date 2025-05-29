from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities and conditions
equations = [
    Eq(D - A, 66),
    Eq(F - G, -145),
    Eq(B - G, -209),
    Eq(E, 3.0 * B),
    Eq(I - A, 112),
    Eq(A + H, 159),
    Eq(D + H, 225),
    Eq(E, 2.4 * J),
    Eq(A - H, -141)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

print(result)