from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, F, G, H, I, J = symbols('A B C D F G H I J')

# Define the equations
eq1 = Eq(C - I, 13)
eq2 = Eq(G - I, 147)
eq3 = Eq(A, 3*D)
eq4 = Eq(G + J, 186)
eq5 = Eq(D + H, 123)
eq6 = Eq(F - A, -216)
eq7 = Eq(H + J, 84)
eq8 = Eq(F + G, 159)
eq9 = Eq(A + B, 346)
eq10 = Eq(F, 3*I)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, F, G, H, I, J))

# Extract the values in alphabetical order, excluding E
result = [solution[A], solution[B], solution[C], solution[D], solution[F], solution[G], solution[H], solution[I], solution[J]]

# Find the remaining number for E
all_numbers = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}
used_numbers = set(result)
E_value = list(all_numbers - used_numbers)[0]

# Append E's value to the result
result.insert(4, E_value)

print(result)