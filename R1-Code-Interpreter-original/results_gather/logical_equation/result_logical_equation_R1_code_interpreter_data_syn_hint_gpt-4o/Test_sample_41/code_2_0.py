from sympy import symbols, Eq, solve

# Define the symbols
B, C, D, F, H, J, K = symbols('B C D F H J K')

# Define the equations based on the inequalities
eq1 = Eq(C, 3 * J)
eq2 = Eq(F + K, 375)
eq3 = Eq(H, 3 * D)
eq4 = Eq(D + F, 234)
eq5 = Eq(J - B, -59)
eq6 = Eq(J + K, 166)
eq7 = Eq(B - H, 48)
eq8 = Eq(K - C, 102)
eq9 = Eq(C + K, 198)
eq10 = Eq(D - B, -66)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (B, C, D, F, H, J, K))

# Extract the values in alphabetical order
# Note: A, E, G, I are not involved in the equations, so we need to assign them manually from the remaining numbers
remaining_numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
used_numbers = [solution[B], solution[C], solution[D], solution[F], solution[H], solution[J], solution[K]]
unused_numbers = [num for num in remaining_numbers if num not in used_numbers]

# Assign unused numbers to A, E, G, I
result = [unused_numbers[0], solution[B], solution[C], solution[D], unused_numbers[1], solution[F], unused_numbers[2], solution[H], unused_numbers[3], solution[J], solution[K]]

print(result)