from sympy import symbols, Eq, solve, And

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations and inequalities
equations = [
    Eq(E + F, 8),
    Eq(C - D, 14),
    Eq(A + M, 35),
    Eq(D + K, 3),
    Eq(G + J, 60),
    Eq(H + I, 60),
    Eq(A - K, 5),
    Eq(G - H, 26),
    Eq(G - F, 45),
    B > F
]

# Define the domain of possible values
possible_values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Solve the system of equations and inequalities
solution = solve(equations + [And(*[var in possible_values for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]])])

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

# Print the result
print(f"<<<{result}>>>")