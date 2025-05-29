from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
eq1 = Eq(D + E, 15)
eq2 = Eq(A + F, 201)
eq3 = Eq(D, 1.5 * E)
eq4 = Eq(C, 4.0 * D)
eq5 = Eq(E - I, -26)
eq6 = Eq(E, 2.0 * G)
eq7 = Eq(C + I, 68)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7), (A, B, C, D, E, F, G, H, I))

# Extract the values from the solution
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]

# Print the values in alphabetical order
print(values)