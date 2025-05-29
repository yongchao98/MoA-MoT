from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
eq1 = Eq(D - K, 23)
eq2 = Eq(L, 3.6 * A)
eq3 = Eq(A + H, 25)
eq4 = Eq(K, 2.5 * C)
eq5 = Eq(E, 2.4 * A)
eq6 = Eq(H, 1.5 * A)
eq7 = Eq(I - G, -42)
eq8 = Eq(I - E, -21)
eq9 = Eq(L, 1.5 * E)

# Solve the equations
solutions = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9), (A, E, H, L, K, C, D, I, G))

# Extract the solution
A_val = solutions[A]
E_val = solutions[E]
H_val = solutions[H]
L_val = solutions[L]
K_val = solutions[K]
C_val = solutions[C]
D_val = solutions[D]
I_val = solutions[I]
G_val = solutions[G]

# Remaining numbers to assign: 1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50
# Assigned numbers: A_val, E_val, H_val, L_val, K_val, C_val, D_val, I_val, G_val

# Find the remaining numbers
all_numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
assigned_numbers = {A_val, E_val, H_val, L_val, K_val, C_val, D_val, I_val, G_val}
remaining_numbers = list(all_numbers - assigned_numbers)

# Assign remaining numbers to B, F, J, M ensuring B > F
B_val, F_val, J_val, M_val = sorted(remaining_numbers, reverse=True)

# Ensure B > F
if B_val < F_val:
    B_val, F_val = F_val, B_val

# Return the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, M_val]
print(result)