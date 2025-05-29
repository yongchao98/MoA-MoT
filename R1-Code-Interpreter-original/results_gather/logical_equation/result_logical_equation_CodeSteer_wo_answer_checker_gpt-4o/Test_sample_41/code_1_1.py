from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given inequalities
equations = [
    Eq(C, 3 * J),        # C = 3.0J
    Eq(F + K, 375),      # F + K = 375
    Eq(H, 3 * D),        # H = 3.0D
    Eq(D + F, 234),      # D + F = 234
    Eq(J - B, -59),      # J - B = -59
    Eq(J + K, 166),      # J + K = 166
    Eq(B - H, 48),       # B - H = 48
    Eq(K - C, 102),      # K - C = 102
    Eq(C + K, 198),      # C + K = 198
    Eq(D - B, -66)       # D - B = -66
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Given numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Map the solution to the given numbers
result = [solution[var] for var in [A, B, C, D, E, F, G, H, I, J, K]]

# Output the results in alphabetical order
print(f"<<<{result}>>>")