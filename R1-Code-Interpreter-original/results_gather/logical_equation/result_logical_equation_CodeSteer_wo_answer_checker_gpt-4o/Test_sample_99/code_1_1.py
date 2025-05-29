from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations
equations = [
    Eq(B + I, 31),
    Eq(I, 1.5 * D),
    Eq(K - E, -2),
    Eq(B + H, 57),
    Eq(C + F, 132),
    Eq(K - D, -13),
    Eq(M, 2.0 * E),
    Eq(A + G, 43),
    Eq(D - H, -34),
    Eq(F, 4.0 * I),
    Eq(I + L, 69)
]

# Define the inequality
inequality = M > B

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations and inequality
solution = solve(equations + [inequality], (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Filter solutions to match possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
        print(f"<<<{result}>>>")
        break