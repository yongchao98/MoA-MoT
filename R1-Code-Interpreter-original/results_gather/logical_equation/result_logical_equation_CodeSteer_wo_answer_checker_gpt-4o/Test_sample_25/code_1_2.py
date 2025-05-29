from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Set up the equations based on the given inequalities
equations = [
    Eq(D, 1.5 * I),
    Eq(I, 2.5 * H),
    Eq(C + H, 52),
    Eq(I, 4.0 * C),
    Eq(C - H, -12),
    Eq(D + I, 200),
    Eq(B - H, -23),
    Eq(H, 1.6 * C)
]

# Solve the system of equations
solution = solve(equations, (C, D, H, I, B))

# Extract the values from the solution
C_value = solution[C]
H_value = solution[H]
I_value = solution[I]
D_value = solution[D]
B_value = solution[B]

# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Used numbers from the solution
used_numbers = {C_value, H_value, I_value, D_value, B_value}

# Assign the remaining numbers to A, E, F, and G
remaining_numbers = sorted(set(numbers) - used_numbers)
A_value, E_value, F_value, G_value = remaining_numbers

# Return the values in alphabetical order
result = [A_value, B_value, C_value, D_value, E_value, F_value, G_value, H_value, I_value]
print(f"<<<{result}>>>")