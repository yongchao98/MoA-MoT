from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(I - E, -4)
eq2 = Eq(F, 2.5 * K)
eq3 = Eq(I, 1.5 * B)
eq4 = Eq(L - A, 33)
eq5 = Eq(E, 2.8 * H)
eq6 = Eq(C + M, 111)
eq7 = Eq(L - C, -60)
eq8 = Eq(H, 2.0 * F)
eq9 = Eq(F + L, 41)
eq10 = Eq(K - E, -26)
eq11 = Eq(A + F, 8)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the values in alphabetical order
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(values)