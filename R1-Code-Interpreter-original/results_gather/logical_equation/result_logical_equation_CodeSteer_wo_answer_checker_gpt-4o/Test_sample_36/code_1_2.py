from sympy import symbols, Eq, nonlinsolve

# Define the variables
A, B, C, E, F, H, I, J, K = symbols('A B C E F H I J K')

# Define the equations
equations = [
    Eq(F, 4 * C),
    Eq(J - A, -48),
    Eq(K - B, -75),
    Eq(B, 3 * A),
    Eq(B + C, 234),
    Eq(F + J, 63),
    Eq(H + K, 271),
    Eq(A - I, 59),
    Eq(C, 3 * E)
]

# Define the possible values
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solutions = nonlinsolve(equations, (A, B, C, E, F, H, I, J, K))

# Filter solutions to match the possible values
for solution in solutions:
    if all(value in possible_values for value in solution):
        # Assign the remaining numbers to D and G
        remaining_values = set(possible_values) - set(solution)
        D, G = remaining_values
        result = [solution[A], solution[B], solution[C], D, solution[E], solution[F], G, solution[H], solution[I], solution[J], solution[K]]
        print(f"<<<{result}>>>")
        break