from sympy import symbols, Eq, solve

# Define symbolic variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(F + H, 6),
    Eq(E, 2.4 * D),
    Eq(B, 2.0 * F),
    Eq(C, 3.0 * M),
    Eq(I - B, 34),
    Eq(H, 2.5 * B),
    Eq(E + M, 39),
    Eq(L, 4.0 * J),
    Eq(M, 3.0 * H),
    Eq(L, 2.8 * D)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K, L, M), dict=True)

# Extract the solution
if solution:
    solution = solution[0]
    # Assign the remaining numbers
    solution[A] = 3
    solution[G] = 16
    solution[K] = 50

    # Create a list of values in alphabetical order
    result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F],
              solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

    # Print the result
    print(f"<<<{result}>>>")
else:
    print("No solution found.")