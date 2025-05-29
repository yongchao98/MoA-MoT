from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the given conditions
equations = [
    Eq(G, 1.5 * H),
    Eq(E - D, -41),
    Eq(F, 1.5 * E),
    Eq(D + G, 130),
    Eq(F + H, 126),
    Eq(F + I, 123),
    Eq(E, 2.5 * B),
    Eq(D - E, 41)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I))

# Check if the solution is a list or a dictionary
if isinstance(solution, list):
    # If it's a list, take the first solution (assuming it's the correct one)
    solution = solution[0]

# Extract the solution and print it in the required format
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I]]
print(f"<<<{result}>>>")