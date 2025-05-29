from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(K, 2.0 * E)
eq2 = Eq(E - C, 72)
eq3 = Eq(J - G, -205)
eq4 = Eq(A + F, 128)
eq5 = Eq(C + F, 51)
eq6 = Eq(G, 1.5 * K)
eq7 = Eq(B + C, 19)
eq8 = Eq(H, 3.0 * C)
eq9 = Eq(F, 2.4 * J)
eq10 = Eq(G + J, 245)
ineq1 = K > H

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, ineq1), (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Extract the solution
solution_values = solution[0]

# List of numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Assign numbers to the letters based on the solution
assigned_values = [solution_values[A], solution_values[B], solution_values[C], solution_values[D], solution_values[E], solution_values[F], solution_values[G], solution_values[H], solution_values[I], solution_values[J], solution_values[K]]

# Print the assigned values
print(assigned_values)