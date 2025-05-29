from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
eq1 = Eq(G, 3.6 * J)
eq2 = Eq(C, 3.0 * E)
eq3 = Eq(E + I, 18)
eq4 = Eq(K - C, 5)
eq5 = Eq(L - B, 17)
eq6 = Eq(B + H, 12)
eq7 = Eq(B + C, 52)
eq8 = Eq(D - E, -13)
eq9 = Eq(B - A, -9)
eq10 = Eq(A - J, 6)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the values in alphabetical order
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(values)