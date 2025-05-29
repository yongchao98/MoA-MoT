from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, G, H, I, J, K = symbols('A B C D E G H I J K')

# Define the equations
eq1 = Eq(B + E, 170)
eq2 = Eq(H - C, -46)
eq3 = Eq(K, 3.0 * G)
eq4 = Eq(B + I, 230)
eq5 = Eq(C - H, 46)
eq6 = Eq(D, 3.0 * H)
eq7 = Eq(C - D, -104)
eq8 = Eq(H - J, 36)
eq9 = Eq(A, 2.4 * E)
eq10 = Eq(I - E, 60)
eq11 = Eq(C - K, 112)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, C, D, E, G, H, I, J, K))

# Output the solution in the specified format
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[G], solution[H], solution[I], solution[J], solution[K]]
print(f"<<<{result}>>>")