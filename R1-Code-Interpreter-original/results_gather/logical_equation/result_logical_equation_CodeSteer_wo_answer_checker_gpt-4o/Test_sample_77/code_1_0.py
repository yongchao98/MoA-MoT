from sympy import symbols, Eq, solve

# Define the variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations and inequalities
equations = [
    Eq(K - G, 25),
    Eq(F, 3.5 * M),
    Eq(B + E, 141),
    Eq(J + K, 38),
    Eq(B + I, 50),
    Eq(A - G, 33),
    Eq(I + M, 7),
    Eq(A - H, -14),
    Eq(M - H, -48),
    Eq(C - E, -72),
    Eq(F, 1.4 * I)
]

# Define the possible values for each variable
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Function to check if a solution matches the possible values
def is_valid_solution(sol):
    return all(sol[var] in possible_values for var in sol)

# Filter solutions to match possible values
valid_solutions = [sol for sol in solution if is_valid_solution(sol)]

# Print the first valid solution
if valid_solutions:
    sol = valid_solutions[0]
    result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K], sol[L], sol[M]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")