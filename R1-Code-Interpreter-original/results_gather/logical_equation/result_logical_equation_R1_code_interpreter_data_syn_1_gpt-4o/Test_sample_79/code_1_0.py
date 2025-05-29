from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations
equations = [
    Eq(F + J, 15),
    Eq(F - G, -91),
    Eq(D, 3.5 * K),
    Eq(E + G, 141),
    Eq(A - K, 13),
    Eq(F + M, 29),
    Eq(D, 1.4 * F),
    Eq(L, 4.0 * D),
    Eq(C, 3.6 * J),
    Eq(A + M, 39),
    Eq(C, 2.4 * A)
]

# Solve the equations
solution = solve(equations, (A, C, D, E, F, G, J, K, L, M))

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Extract the solution in alphabetical order for the solved variables
solved_values = [solution[A], solution[C], solution[D], solution[E], solution[F], solution[G], solution[J], solution[K], solution[L], solution[M]]

# Determine the remaining letters (B, H, I) by elimination
remaining_letters = ['B', 'H', 'I']
remaining_values = [num for num in numbers if num not in solved_values]

# Combine the results
result = [solution[A], remaining_values[0], solution[C], solution[D], solution[E], solution[F], solution[G], remaining_values[1], remaining_values[2], solution[J], solution[K], solution[L], solution[M]]

print(result)