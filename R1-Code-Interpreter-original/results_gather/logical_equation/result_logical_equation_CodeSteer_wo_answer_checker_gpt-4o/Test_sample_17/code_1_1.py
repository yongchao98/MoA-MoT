from sympy import symbols, Eq, solve

# Define the symbols for each letter
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
equations = [
    Eq(A, 1.6 * G),
    Eq(E + H, 129),
    Eq(H, 1.5 * B),
    Eq(E + I, 130),
    Eq(A + B, 112),
    Eq(D + G, 23),
    Eq(B, 2.5 * A),
    Eq(B, 4.0 * G)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], 36, solution[D], solution[E], 6, solution[G], solution[H], solution[I]]

# Print the result in the required format
print(f"<<<{result}>>>")