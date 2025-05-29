from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(C - G, -36)
eq2 = Eq(C + G, 42)
eq3 = Eq(C - F, -17)
eq4 = Eq(F - G, -19)
eq5 = Eq(B + F, 68)
eq6 = Eq(E + H, 137)
eq7 = Eq(C + H, 124)
eq8 = Eq(A + K, 375)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8), (A, B, C, E, F, G, H, K))

# Extract the values
C_value = solution[C]
G_value = solution[G]
F_value = solution[F]
B_value = solution[B]
E_value = solution[E]
H_value = solution[H]
A_value = solution[A]
K_value = solution[K]

# Find the remaining values
remaining_numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
used_numbers = [C_value, G_value, F_value, B_value, E_value, H_value, A_value, K_value]
I_value = next(num for num in remaining_numbers if num not in used_numbers)
J_value = next(num for num in remaining_numbers if num not in used_numbers and num > C_value)

# Print the results in alphabetical order
print([A_value, B_value, C_value, D, E_value, F_value, G_value, H_value, I_value, J_value, K_value])