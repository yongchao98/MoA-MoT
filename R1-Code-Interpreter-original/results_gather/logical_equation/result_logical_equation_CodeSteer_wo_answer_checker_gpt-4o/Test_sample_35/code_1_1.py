from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the constraints
equations = [
    Eq(E, 1.5 * K),
    Eq(D + I, 155),
    Eq(K, 2.0 * D),
    Eq(E, 3.0 * D),
    Eq(D - C, -46),
    Eq(F + H, 51),
    Eq(H + I, 83),
    Eq(F, 3.0 * G),
    Eq(E + K, 375),
    Eq(F - D, -27)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution values
# The solution will be a dictionary with variable keys and their corresponding values
values = [solution[var] for var in [A, B, C, D, E, F, G, H, I, J, K]]

# Output the values in alphabetical order
print(f"<<<{values}>>>")