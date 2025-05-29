from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations and inequalities
equations = [
    Eq(D - A, 66),
    Eq(F - G, -145),
    Eq(B - G, -209),
    Eq(E, 3.0 * B),
    Eq(I - A, 112),
    Eq(A + H, 159),
    Eq(D + H, 225),
    Eq(E, 2.4 * J),
    Eq(A - H, -141)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter the solution to match the possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
        break

# Print the result
print(f"<<<{result}>>>")