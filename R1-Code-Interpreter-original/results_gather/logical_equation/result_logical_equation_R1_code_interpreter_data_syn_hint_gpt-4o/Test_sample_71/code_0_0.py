from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
equations = [
    Eq(B - F, -21),
    Eq(F + L, 73),
    Eq(F - A, -22),
    Eq(K - L, -42),
    Eq(B + C, 43),
    Eq(H + L, 60),
    Eq(H + M, 17),
    Eq(C, 2.4 * H),
    Eq(D + K, 4),
    Eq(G, 2.5 * M)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(result)