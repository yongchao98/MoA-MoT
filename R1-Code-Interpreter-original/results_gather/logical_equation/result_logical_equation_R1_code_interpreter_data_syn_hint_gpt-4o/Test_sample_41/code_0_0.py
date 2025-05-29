from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities
eq1 = Eq(C, 3 * J)
eq2 = Eq(F + K, 375)
eq3 = Eq(H, 3 * D)
eq4 = Eq(D + F, 234)
eq5 = Eq(J - B, -59)
eq6 = Eq(J + K, 166)
eq7 = Eq(B - H, 48)
eq8 = Eq(K - C, 102)
eq9 = Eq(C + K, 198)
eq10 = Eq(D - B, -66)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

print(result)