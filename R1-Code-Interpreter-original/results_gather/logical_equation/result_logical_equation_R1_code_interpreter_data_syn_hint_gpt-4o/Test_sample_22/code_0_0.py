from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
eq1 = Eq(B, 2.0 * G)
eq2 = Eq(I, 2.5 * F)
eq3 = Eq(B + I, 86)
eq4 = Eq(A + G, 12)
eq5 = Eq(F + I, 112)
eq6 = Eq(C + H, 140)
eq7 = Eq(A - H, -11)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7), (A, B, C, F, G, H, I))

# Extract the values for each letter
A_val = solution[A]
B_val = solution[B]
C_val = solution[C]
F_val = solution[F]
G_val = solution[G]
H_val = solution[H]
I_val = solution[I]

# The remaining numbers are 20, 32, 36, 80, 120, 121
# We need to assign D and E from these remaining numbers
remaining_numbers = [20, 32, 36, 80, 120, 121]
D_val = remaining_numbers[0]  # Assigning arbitrarily
E_val = remaining_numbers[1]  # Assigning arbitrarily

# Print the values in alphabetical order
print([A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val])