from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(J, 3.5 * D),
    Eq(K - F, -12),
    Eq(G + M, 26),
    Eq(J + K, 31),
    Eq(C - D, 13),
    Eq(G, 1.6 * M),
    Eq(A + L, 124),
    Eq(M - A, -18),
    Eq(C, 3.0 * H),
    Eq(J > D),
    Eq(L > G)
]

# Possible values for the letters
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter solutions to match the possible values
for sol in solution:
    if all(sol[var] in values for var in sol):
        print([sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]])
        break