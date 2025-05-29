from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities
equations = [
    Eq(C + J, 186),
    Eq(H, 3 * F),
    Eq(C - B, -44),
    Eq(E + H, 234),
    Eq(E + I, 130),
    Eq(I - J, -29),
    Eq(D + E, 57),
    Eq(K - B, -64),
    Eq(C + D, 84),
    Eq(B - G, 53)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

print(result)