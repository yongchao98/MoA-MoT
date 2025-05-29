from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities
eq1 = Eq(K - H, -150)
eq2 = Eq(G - K, 46)
eq3 = Eq(B - F, 70)
eq4 = Eq(K - A, 59)
eq5 = Eq(H, 3.0 * K)
eq6 = Eq(C + H, 252)
eq7 = Eq(D, 4.0 * I)
eq8 = Eq(C + F, 107)
eq9 = Eq(J - B, -102)
eq10 = Eq(G - I, 112)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

print(result)