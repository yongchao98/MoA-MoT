from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities
eq1 = Eq(K - H, -150)
eq2 = Eq(G - K, 46)
eq3 = Eq(B - F, 70)
eq4 = Eq(K - A, 59)
eq5 = Eq(H, 3 * K)
eq6 = Eq(C + H, 252)
eq7 = Eq(D, 4 * I)
eq8 = Eq(C + F, 107)
eq9 = Eq(J - B, -102)
eq10 = Eq(G - I, 112)

# List of possible values
values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, E, F, G, H, I, J, K))

# Filter the solution to match the possible values
final_solution = {var: val for var, val in solution.items() if val in values}

# Sort the solution by the variable names
sorted_solution = [final_solution[var] for var in sorted(final_solution)]

print(sorted_solution)