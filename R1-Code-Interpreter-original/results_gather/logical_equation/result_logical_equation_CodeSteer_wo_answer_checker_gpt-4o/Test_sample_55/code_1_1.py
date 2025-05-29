from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given constraints
equations = [
    Eq(C - K, -17),
    Eq(J, 2 * D),
    Eq(G - K, 205),
    Eq(E + H, 119),
    Eq(C - J, -147),
    Eq(F + K, 36),
    Eq(F + G, 241),
    Eq(H - F, 64),
    Eq(I - B, -39),
    Eq(D + K, 95),
    Eq(B + J, 198)
]

# Define the possible values for each variable
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter the solution to match the possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
        break

# Print the result in the required format
print(f"<<<{result}>>>")