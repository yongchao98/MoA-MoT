from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the problem statement
equations = [
    Eq(D + H, 13),
    Eq(E - C, 15),
    Eq(C + H, 5),
    Eq(H, 4 * C),
    Eq(H - B, -117),
    Eq(F, 3 * C),
    Eq(F + H, 7)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I))

# Extract the solution
# Since the solution is a dictionary, we can extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]

# Print the result
print(f"<<<{result}>>>")