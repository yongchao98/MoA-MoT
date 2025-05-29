from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(B - F, -21)
eq2 = Eq(F + L, 73)
eq3 = Eq(F - A, -22)
eq4 = Eq(K - L, -42)
eq5 = Eq(B + C, 43)
eq6 = Eq(H + L, 60)
eq7 = Eq(H + M, 17)
eq8 = Eq(C, 2.4 * H)
eq9 = Eq(D + K, 4)
eq10 = Eq(G, 2.5 * M)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(result)