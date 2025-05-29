from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the inequalities
equations = [
    Eq(E - D, 25),
    Eq(H + I, 103),
    Eq(E, 2.8 * M),
    Eq(C + I, 120),
    Eq(A + I, 101),
    Eq(I - J, 51),
    Eq(I + J, 141),
    Eq(H - I, -89),
    Eq(E, 4.0 * H),
    Eq(G, 3.0 * A),
    Eq(G - E, -13)
]

# Solve the equations for the known variables
solution = solve(equations, (A, C, D, E, G, H, I, J, M))

# Assign the remaining numbers to the remaining variables
# Known numbers: 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96
# Assign the remaining numbers to B, F, K, L
remaining_numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
assigned_numbers = [solution[A], solution[C], solution[D], solution[E], solution[G], solution[H], solution[I], solution[J], solution[M]]
remaining_numbers = [num for num in remaining_numbers if num not in assigned_numbers]

# Assign the remaining numbers to B, F, K, L
solution[B] = remaining_numbers[0]
solution[F] = remaining_numbers[1]
solution[K] = remaining_numbers[2]
solution[L] = remaining_numbers[3]

# Extract the solution and sort it in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(f"<<<{result}>>>")