from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(A, 2.8 * K)
eq2 = Eq(A + F, 30)
eq3 = Eq(C + D, 31)
eq4 = Eq(G - L, 45)
eq5 = Eq(L - C, -19)
eq6 = Eq(I, 1.5 * K)
eq7 = Eq(J + K, 13)
eq8 = Eq(B, 3.6 * K)
eq9 = Eq(C, 1.6 * I)
eq10 = Eq(B + G, 86)
eq11 = Eq(K - J, 7)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, C, D, F, G, I, J, K, L))

# Extract the values
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

# Print the values
print(values)