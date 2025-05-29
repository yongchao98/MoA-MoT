from sympy import symbols, Eq, solve

# Define symbolic variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
equations = [
    Eq(A, 3 * E),      # A = 3.0E
    Eq(D + F, 81),     # D + F = 81
    Eq(A + D, 28),     # A + D = 28
    Eq(C + F, 96),     # C + F = 96
    Eq(A + E, 36),     # A + E = 36
    Eq(B + E, 12),     # B + E = 12
    Eq(B - F, -77)     # B - F = -77
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F))

# Extract the values for A, B, C, D, E, F
A_val = solution[A]
B_val = solution[B]
C_val = solution[C]
D_val = solution[D]
E_val = solution[E]
F_val = solution[F]

# Remaining numbers to be assigned to G, H, I
remaining_numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
used_numbers = [A_val, B_val, C_val, D_val, E_val, F_val]
remaining_numbers = [num for num in remaining_numbers if num not in used_numbers]

# Assign the remaining numbers to G, H, I
G_val, H_val, I_val = remaining_numbers

# Output the values in alphabetical order
result = [A_val, B_val, C_val, D_val, E_val, F_val, G_val, H_val, I_val]
print(f"<<<{result}>>>")