from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(D, A + 43),
    Eq(F, 1.5 * C),
    Eq(A + F, 43),
    Eq(C, 1.6 * L),
    Eq(L, 3.0 * J),
    Eq(F - J, 31),
    Eq(M, 4.5 * G),
    Eq(E + G, 12),
    Eq(J + K, 33),
    Eq(F, 3.6 * G)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(result)