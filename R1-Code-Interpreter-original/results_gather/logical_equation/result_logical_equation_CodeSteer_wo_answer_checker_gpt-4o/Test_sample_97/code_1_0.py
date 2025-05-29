from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given inequalities
equations = [
    Eq(A - L, -5),
    Eq(I, 3.2 * C),
    Eq(J + M, 60),
    Eq(L + M, 31),
    Eq(C, 2.5 * A),
    Eq(L, 3.5 * A),
    Eq(D, 4.0 * M),
    Eq(B, 4.5 * G),
    Eq(I + M, 40),
    Eq(A + D, 98),
    Eq(F, 1.5 * G),
    Eq(D + G, 106)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, F, G, I, J, L, M))

# Extract the values from the solution
A_val = solution[A]
B_val = solution[B]
C_val = solution[C]
D_val = solution[D]
F_val = solution[F]
G_val = solution[G]
I_val = solution[I]
J_val = solution[J]
L_val = solution[L]
M_val = solution[M]

# Remaining numbers to assign: 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96
# Assign remaining numbers to E, H, and K
remaining_numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
assigned_numbers = [A_val, B_val, C_val, D_val, F_val, G_val, I_val, J_val, L_val, M_val]
unassigned_numbers = list(set(remaining_numbers) - set(assigned_numbers))

E_val = unassigned_numbers[0]
H_val = unassigned_numbers[1]
K_val = unassigned_numbers[2]

# Output the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val, J_val, K_val, L_val, M_val]
print(f"<<<{result}>>>")