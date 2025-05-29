from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the inequalities and equalities
equations = [
    Eq(H, 2.4 * A),
    Eq(H, 3.0 * F),
    Eq(B - I, -216),
    Eq(C - J, -72),
    Eq(I, 1.5 * G),
    Eq(I, 3.0 * J),
    Eq(B - H, -39),
    Eq(D - I, -104),
    Eq(J - F, 59),
    Eq(B - F, -7)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution and sort it in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

print(result)