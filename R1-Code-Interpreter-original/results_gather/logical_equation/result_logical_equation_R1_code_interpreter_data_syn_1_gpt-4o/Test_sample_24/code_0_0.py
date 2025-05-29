from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
eq1 = Eq(G, 1.5 * H)
eq2 = Eq(E - D, -41)
eq3 = Eq(F, 1.5 * E)
eq4 = Eq(D + G, 130)
eq5 = Eq(F + H, 126)
eq6 = Eq(F + I, 123)
eq7 = Eq(E, 2.5 * B)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7), (A, B, C, D, E, F, G, H, I))

# Extract the solution and sort it in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]
print(result)