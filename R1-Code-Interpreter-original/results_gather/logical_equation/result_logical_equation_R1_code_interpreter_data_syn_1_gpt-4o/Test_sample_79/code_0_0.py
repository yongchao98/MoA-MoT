from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
equations = [
    Eq(F + J, 15),
    Eq(F - G, -91),
    Eq(D, 3.5 * K),
    Eq(E + G, 141),
    Eq(A - K, 13),
    Eq(F + M, 29),
    Eq(D, 1.4 * F),
    Eq(L, 4.0 * D),
    Eq(C, 3.6 * J),
    Eq(A + M, 39),
    Eq(C, 2.4 * A)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(result)