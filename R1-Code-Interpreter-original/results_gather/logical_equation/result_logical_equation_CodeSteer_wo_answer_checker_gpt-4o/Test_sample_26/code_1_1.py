from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the constraints
equations = [
    Eq(B - I, -74),
    Eq(C + E, 140),
    Eq(A + G, 39),
    Eq(E + H, 29),
    Eq(B, 2 * A),
    Eq(E + I, 100),
    Eq(A + F, 124),
    Eq(A - H, -6)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I))

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]

# Print the result in the required format
print(f"<<<{result}>>>")