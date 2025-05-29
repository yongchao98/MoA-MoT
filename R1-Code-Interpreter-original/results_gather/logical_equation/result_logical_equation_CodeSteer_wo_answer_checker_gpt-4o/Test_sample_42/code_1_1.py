from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given inequalities
equations = [
    Eq(H + J, 130),
    Eq(C + I, 177),
    Eq(C + G, 75),
    Eq(G, 3 * F),
    Eq(B - E, 150),
    Eq(J - C, -18),
    Eq(E + G, 123),
    Eq(D + F, 96),
    Eq(E + K, 78),
    Eq(H + K, 124)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution and print it in the required format
values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]
print(f"<<<{values}>>>")