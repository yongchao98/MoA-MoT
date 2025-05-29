from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(E + K, 81)
eq2 = Eq(E, 1.5 * G)
eq3 = Eq(B + J, 5)
eq4 = Eq(M - F, -6)
eq5 = Eq(K, 4.5 * M)
eq6 = Eq(F, 3.2 * C)
eq7 = Eq(C - E, -31)
eq8 = Eq(L, 1.4 * C)
eq9 = Eq(A + M, 38)
eq10 = Eq(F + I, 112)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, E, F, G, I, K, L, M))

# Extract the solution
solution_values = [solution[A], solution[B], solution[C], solution[E], solution[F], solution[G], solution[I], solution[K], solution[L], solution[M]]

# Print the solution
print(solution_values)