from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(F - C, -141)
eq2 = C > D
eq3 = Eq(E + K, 305)
eq4 = Eq(K, 3.0 * A)
eq5 = Eq(H, 3.0 * F)
eq6 = Eq(I + J, 157)
eq7 = Eq(K, 1.5 * C)
eq8 = Eq(K - A, 150)
eq9 = Eq(F + G, 25)
eq10 = Eq(B + C, 153)

# List of possible values
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve((eq1, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter solutions that satisfy all conditions
valid_solution = None
for sol in solution:
    if all(val in values for val in sol.values()) and sol[C] > sol[D]:
        valid_solution = sol
        break

# Extract the values in alphabetical order
result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D], valid_solution[E],
          valid_solution[F], valid_solution[G], valid_solution[H], valid_solution[I], valid_solution[J], valid_solution[K]]

print(result)