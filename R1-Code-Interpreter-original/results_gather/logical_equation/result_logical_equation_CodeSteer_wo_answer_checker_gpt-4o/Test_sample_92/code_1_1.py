from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations and inequalities
equations = [
    Eq(A, 3.0 * L),
    Eq(A + M, 39),
    Eq(L + M, 29),
    Eq(B, 3.2 * L),
    Eq(A, 1.5 * I),
    Eq(C - B, 80),
    Eq(E + K, 43),
    C > D,
    Eq(E, 2.4 * A),
    Eq(C + D, 124),
    Eq(H - A, -12),
    Eq(F + G, 47)
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter the solution to match the possible values
for sol in solution:
    if all(value in possible_values for value in sol.values()):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        print(f"<<<{result}>>>")
        break