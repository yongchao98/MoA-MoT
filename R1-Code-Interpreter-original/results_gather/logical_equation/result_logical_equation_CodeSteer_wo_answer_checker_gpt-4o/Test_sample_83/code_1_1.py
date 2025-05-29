from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the inequalities
equations = [
    Eq(E - D, 25),
    Eq(H + I, 103),
    Eq(E, 2.8 * M),
    Eq(C + I, 120),
    Eq(A + I, 101),
    Eq(I - J, 51),
    Eq(I + J, 141),
    Eq(H - I, -89),
    Eq(E, 4.0 * H),
    Eq(G, 3.0 * A),
    Eq(G - E, -13)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M))

# Extract the solution and sort it in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(f"<<<{result}>>>")