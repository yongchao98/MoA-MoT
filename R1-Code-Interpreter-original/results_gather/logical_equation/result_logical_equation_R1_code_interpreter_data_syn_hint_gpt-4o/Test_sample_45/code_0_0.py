from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(F, 2.4 * D)
eq2 = Eq(G, 1.5 * K)
eq3 = Eq(G + K, 375)
eq4 = Eq(F - G, -177)
eq5 = Eq(A + D, 29)
eq6 = Eq(I - D, 19)
eq7 = Eq(A, 3.0 * H)
eq8 = Eq(F, 3.0 * B)
eq9 = Eq(F + K, 198)
eq11 = Eq(J - F, 73)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq11), (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

# Print the values
print(values)