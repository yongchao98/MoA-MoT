from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the problem statement
equations = [
    Eq(K - H, -150),
    Eq(G - K, 46),
    Eq(B - F, 70),
    Eq(K - A, 59),
    Eq(H, 3 * K),
    Eq(C + H, 252),
    Eq(D, 4 * I),
    Eq(C + F, 107),
    Eq(J - B, -102),
    Eq(G - I, 112)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

# Print the result in the required format
print(f"<<<{result}>>>")