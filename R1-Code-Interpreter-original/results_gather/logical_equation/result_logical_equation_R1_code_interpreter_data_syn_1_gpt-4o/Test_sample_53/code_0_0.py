from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(F - K, 36)
eq2 = Eq(B - K, 111)
eq3 = Eq(B + E, 230)
eq4 = Eq(H + J, 25)
eq5 = Eq(E - A, 77)
eq6 = Eq(H + I, 36)
eq7 = Eq(C + E, 305)
eq8 = Eq(C - I, 205)
eq9 = Eq(C + J, 234)
eq10 = Eq(G, 3 * H)
eq11 = Eq(A + K, 42)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, C, D, E, F, G, H, I, J, K))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]
print(result)