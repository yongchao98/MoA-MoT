from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, F, G, I, J, L, M = symbols('A B C D F G I J L M')

# Define the equations based on the constraints
eq1 = Eq(A - L, -5)
eq2 = Eq(I, 3.2 * C)
eq3 = Eq(J + M, 60)
eq4 = Eq(L + M, 31)
eq5 = Eq(C, 2.5 * A)
eq6 = Eq(L, 3.5 * A)
eq7 = Eq(D, 4.0 * M)
eq8 = Eq(B, 4.5 * G)
eq9 = Eq(I + M, 40)
eq10 = Eq(A + D, 98)
eq11 = Eq(F, 1.5 * G)
eq12 = Eq(D + G, 106)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12), (A, B, C, D, F, G, I, J, L, M))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], None, solution[F], solution[G], None, solution[I], solution[J], None, solution[L], solution[M]]

print(result)