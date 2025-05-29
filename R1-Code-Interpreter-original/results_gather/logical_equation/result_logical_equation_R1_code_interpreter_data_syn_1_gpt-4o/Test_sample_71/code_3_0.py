from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, F, G, H, K, L, M = symbols('A B C D F G H K L M')

# Define the equations
eq1 = Eq(B - F, -21)
eq2 = Eq(F + L, 73)
eq3 = Eq(F - A, -22)
eq4 = Eq(K - L, -42)
eq5 = Eq(B + C, 43)
eq6 = Eq(H + L, 60)
eq7 = Eq(H + M, 17)
eq8 = Eq(C, 2.4 * H)
eq9 = Eq(D + K, 4)
eq10 = Eq(G, 2.5 * M)

# Solve the system of equations
solution = solve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10), (A, B, C, D, F, G, H, K, L, M))

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], None, solution[F], solution[G], solution[H], None, None, solution[K], solution[L], solution[M]]

# Find the remaining numbers for E, I, and J
used_numbers = set(result)
remaining_numbers = [num for num in numbers if num not in used_numbers]

# Assign remaining numbers to E, I, and J
result[4] = remaining_numbers[0]  # E
result[8] = remaining_numbers[1]  # I
result[9] = remaining_numbers[2]  # J

# Ensure all numbers are unique and from the given set
result = [int(x) if x is not None else None for x in result]

print(result)