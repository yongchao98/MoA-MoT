from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(C + H, 52),
    Eq(G + I, 13),
    Eq(H, 2.4 * E),
    Eq(F + J, 57),
    Eq(A - I, 25),
    Eq(J - I, 4),
    Eq(D + G, 34),
    Eq(F + L, 95),
    Eq(A, 4.0 * J),
    Eq(F + G, 60),
    Eq(J, 1.4 * K)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Possible values for the letters
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Check which solution matches the possible values
for sol in solution:
    values = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
    if sorted(values) == sorted(possible_values):
        print(f"<<<{values}>>>")
        break