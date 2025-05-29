from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
eq1 = Eq(E, 3.6 * B)
eq2 = Eq(G, 4.8 * M)
eq3 = Eq(L - C, 12)
eq4 = Eq(D + H, 66)
eq5 = Eq(L - D, -35)
eq6 = Eq(B - D, -40)
eq7 = Eq(H - G, -8)
eq8 = Eq(C, 1.5 * K)
eq9 = Eq(H, 1.6 * B)
eq10 = Eq(J - E, 9)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (B, C, D, E, G, H, J, K, L, M))

# Extract the values
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

# Print the values
print(values)