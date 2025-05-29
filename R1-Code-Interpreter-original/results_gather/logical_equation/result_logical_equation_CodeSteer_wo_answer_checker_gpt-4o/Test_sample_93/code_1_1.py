from sympy import symbols, Eq, solve, S

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations
equations = [
    Eq(C + E, 141),
    Eq(J + L, 33),
    Eq(K, 2.4 * H),
    Eq(B - J, 10),
    Eq(E - C, 51),
    Eq(E + K, 120),
    Eq(J - H, -5),
    Eq(K - B, 9),
    Eq(B - M, -1),
    Eq(M, 3.2 * J)
]

# Define the inequalities
inequalities = [
    C > F,
    C > M
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations and inequalities
solution = solve(equations + inequalities, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True, domain=S.Integers)

# Filter the solution to match the possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        break

# Print the result in the required format
print(f"<<<{result}>>>")