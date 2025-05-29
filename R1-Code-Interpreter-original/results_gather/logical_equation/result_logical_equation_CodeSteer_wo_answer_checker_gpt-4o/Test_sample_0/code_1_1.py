from sympy import symbols, Eq, solve, S

# Define the symbolic variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I', integer=True)

# Define the equations and inequalities
equations = [
    Eq(A, C - 44),
    Eq(B, C - 79),
    Eq(D, 3 * B),
    Eq(I, 3 * D),
    Eq(A, 4 * I),
    Eq(B, E + 3),
    C > H
]

# Possible values for the letters
possible_values = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Solve the system of equations and inequalities
solution = solve(equations, (A, B, C, D, E, F, G, H, I), dict=True, domain=S.Integers)

# Filter solutions to match the possible values
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I]]
        print(f"<<<{result}>>>")
        break