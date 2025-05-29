from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the problem statement
equations = [
    Eq(A + I, 300),
    Eq(E, 2.4 * C),
    Eq(G, 3.0 * K),
    J > A,
    Eq(A - G, 66),
    Eq(J, 4.0 * C),
    Eq(I, 3.0 * A),
    Eq(F, 2.0 * A),
    Eq(I, 1.5 * F),
    Eq(A - C, 55),
    Eq(H - A, -59)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter the solution to match the possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
        print(f"<<<{result}>>>")
        break