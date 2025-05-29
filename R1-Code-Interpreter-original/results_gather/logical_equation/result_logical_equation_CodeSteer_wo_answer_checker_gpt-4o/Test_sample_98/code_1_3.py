from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the constraints
equations = [
    Eq(A, 2.8 * K),
    Eq(A + F, 30),
    Eq(C + D, 31),
    Eq(G - L, 45),
    Eq(L - C, -19),
    Eq(I, 1.5 * K),
    Eq(J + K, 13),
    Eq(B, 3.6 * K),
    Eq(C, 1.6 * I),
    Eq(B + G, 86),
    Eq(K - J, 7)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, F, G, I, J, K, L))

# Assign the remaining numbers to E, H, M
remaining_numbers = set([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]) - set(solution.values())
solution[E], solution[H], solution[M] = sorted(remaining_numbers)

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F],
          solution[G], solution[H], solution[I], solution[J], solution[K], solution[L], solution[M]]

print(f"<<<{result}>>>")