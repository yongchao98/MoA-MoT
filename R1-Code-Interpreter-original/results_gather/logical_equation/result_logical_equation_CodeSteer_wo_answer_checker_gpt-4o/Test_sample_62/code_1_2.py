from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(G, 3.6 * J),
    Eq(C, 3.0 * E),
    Eq(E + I, 18),
    Eq(K - C, 5),
    Eq(L - B, 17),
    Eq(B + H, 12),
    Eq(B + C, 52),
    Eq(D - E, -13),
    Eq(B - A, -9),
    Eq(A - J, 6)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Given numbers
numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}

# Check if the solution matches the given numbers
for sol in solution:
    if all(value in numbers for value in sol.values()):
        # Print the solution in the required format
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        print(f"<<<{result}>>>")
        break