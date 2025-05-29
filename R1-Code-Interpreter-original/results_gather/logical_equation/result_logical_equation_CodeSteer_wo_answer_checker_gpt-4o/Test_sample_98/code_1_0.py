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

# Extract the solved values
solved_values = {var: solution[var] for var in solution}

# List of all possible numbers
possible_numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Assign the remaining numbers to E, H, M
remaining_numbers = set(possible_numbers) - set(solved_values.values())

# Assign remaining numbers to E, H, M
for num in remaining_numbers:
    if num not in solved_values.values():
        if E not in solved_values:
            solved_values[E] = num
        elif H not in solved_values:
            solved_values[H] = num
        elif M not in solved_values:
            solved_values[M] = num

# Extract the solution in alphabetical order
result = [solved_values[A], solved_values[B], solved_values[C], solved_values[D], solved_values[E], solved_values[F],
          solved_values[G], solved_values[H], solved_values[I], solved_values[J], solved_values[K], solved_values[L], solved_values[M]]

print(f"<<<{result}>>>")