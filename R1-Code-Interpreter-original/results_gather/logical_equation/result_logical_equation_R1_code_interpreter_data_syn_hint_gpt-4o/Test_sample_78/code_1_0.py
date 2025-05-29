from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(H, 2.4 * D),
    Eq(G - A, -9),
    Eq(E, 4.0 * G),
    Eq(F - G, 8),
    Eq(B + E, 78),
    Eq(M - K, 40),
    Eq(D + M, 55),
    Eq(A, 3.2 * K),
    Eq(I - C, 94),
    Eq(G + L, 43),
    Eq(L - A, 20)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Filter the solution to match the given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
valid_solution = [solution[letter] for letter in (A, B, C, D, E, F, G, H, I, J, K, L, M) if solution[letter] in numbers]

print(valid_solution)