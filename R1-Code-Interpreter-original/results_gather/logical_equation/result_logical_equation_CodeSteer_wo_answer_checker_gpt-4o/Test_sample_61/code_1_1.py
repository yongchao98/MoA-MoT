from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(C + M, 43),
    Eq(C, 1.5 * L),
    Eq(L, 2.0 * J),
    Eq(A + K, 81),
    Eq(C - B, -35),
    Eq(H - K, -20),
    Eq(F, 4.8 * J),
    Eq(K - M, 8),
    Eq(D + E, 4),
    Eq(K, 2.4 * C)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Possible values for the letters
possible_values = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}

# Filter the solution to ensure all values are in the possible values set
for sol in solution:
    if all(value in possible_values for value in sol.values()):
        # Print the solution in the required format
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        print(f"<<<{result}>>>")
        break