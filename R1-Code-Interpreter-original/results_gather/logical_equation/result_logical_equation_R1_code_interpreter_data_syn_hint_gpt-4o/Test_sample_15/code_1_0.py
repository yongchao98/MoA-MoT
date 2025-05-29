from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations
eq1 = Eq(C + I, 241)
eq2 = Eq(E - H, -44)
eq3 = Eq(I, 1.5 * H)
eq4 = Eq(D + E, 42)
eq5 = Eq(D + H, 86)
eq6 = Eq(F + G, 35)
eq7 = Eq(D, 2.0 * F)
eq8 = Eq(A + I, 140)

# Solve the system of equations
solution = solve((eq3, eq5, eq7, eq4, eq2, eq1, eq8, eq6), (A, C, D, E, F, G, H, I))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]

print(result)