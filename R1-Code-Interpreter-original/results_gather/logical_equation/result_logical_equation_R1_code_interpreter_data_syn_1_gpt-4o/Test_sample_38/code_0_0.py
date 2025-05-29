from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(C + H, 124)
eq2 = Eq(B + F, 305)
eq3 = Eq(G, 3.0 * J)
eq4 = Eq(H - G, -45)
eq5 = Eq(C + J, 137)
eq6 = Eq(J - A, -11)
eq7 = Eq(E - J, 59)
eq8 = Eq(F, 3.0 * E)
eq9 = Eq(A + H, 30)
eq10 = Eq(A - C, -94)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]
print(result)