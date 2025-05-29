from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
eq1 = Eq(H + K, 273)
eq2 = Eq(B - F, 60)
eq3 = Eq(K - B, 145)
eq4 = Eq(F - I, 11)
eq5 = Eq(D - E, -59)
eq6 = Eq(B, 4.0 * F)
eq7 = Eq(G + H, 51)
eq8 = Eq(K, 3.0 * E)
eq9 = Eq(C + K, 264)
eq10 = Eq(I, 3.0 * G)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (B, C, D, E, F, G, H, I, K))

# Extract the values
B_val = solution[B]
C_val = solution[C]
D_val = solution[D]
E_val = solution[E]
F_val = solution[F]
G_val = solution[G]
H_val = solution[H]
I_val = solution[I]
K_val = solution[K]

# Remaining numbers
remaining_numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
used_numbers = [B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, K_val]
A_val = list(set(remaining_numbers) - set(used_numbers))[0]
J_val = list(set(remaining_numbers) - set(used_numbers) - {A_val})[0]

# Return the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val]
print(result)