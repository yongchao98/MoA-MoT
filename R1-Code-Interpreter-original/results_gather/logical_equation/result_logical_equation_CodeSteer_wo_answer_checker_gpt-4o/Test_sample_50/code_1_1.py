from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given conditions
equations = [
    Eq(B + F, 346),
    Eq(A + E, 155),
    Eq(F + G, 228),
    Eq(A + H, 128),
    Eq(A, 4 * D),
    Eq(I - K, -23),
    Eq(J, 2 * E),
    Eq(F, 3 * E),
    Eq(B - J, -29),
    Eq(A + K, 119),
    Eq(H, 3 * I)
]

# Given numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution and map it to the given numbers
result = [solution[var] for var in [A, B, C, D, E, F, G, H, I, J, K]]

# Ensure the solution uses all the given numbers
if sorted(result) == sorted(numbers):
    print(f"<<<{result}>>>")
else:
    print("No valid solution found that uses all the given numbers.")