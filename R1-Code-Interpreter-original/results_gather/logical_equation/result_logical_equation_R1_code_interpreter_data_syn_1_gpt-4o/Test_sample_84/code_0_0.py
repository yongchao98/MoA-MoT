from sympy import symbols, Eq, solve

# Define the variables
A, B, G, J, K, F, M, I, L = symbols('A B G J K F M I L')

# Define the equations
eq1 = Eq(B, 1.4 * G)
eq2 = Eq(A, 3.2 * G)
eq3 = Eq(J, 3.0 * G)
eq4 = Eq(A + J, 31)
eq5 = Eq(A + B, 23)
eq6 = Eq(F, 4.8 * G)
eq7 = Eq(A, 1.6 * K)
eq8 = Eq(F, 2.4 * K)
eq9 = Eq(M - B, 89)
eq10 = Eq(I - B, 43)
eq11 = Eq(B, 3.5 * L)

# Solve the equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11), (A, B, G, J, K, F, M, I, L))

# Extract the values
A_val = solution[A]
B_val = solution[B]
G_val = solution[G]
J_val = solution[J]
K_val = solution[K]
F_val = solution[F]
M_val = solution[M]
I_val = solution[I]
L_val = solution[L]

# Remaining numbers to assign: C, D, E, H
remaining_numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
assigned_numbers = [A_val, B_val, G_val, J_val, K_val, F_val, M_val, I_val, L_val]
remaining_numbers = [num for num in remaining_numbers if num not in assigned_numbers]

# Assign remaining numbers to C, D, E, H
C_val, D_val, E_val, H_val = remaining_numbers[:4]

# Return the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, M_val]
print(result)