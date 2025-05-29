from sympy import symbols, Eq, solve

# Define the symbols
A, B, C, D, F, G, H, I, J, K = symbols('A B C D F G H I J K')

# Define the equations based on the inequalities and equalities
equations = [
    Eq(H, 2.4 * A),
    Eq(H, 3.0 * F),
    Eq(B - I, -216),
    Eq(C - J, -72),
    Eq(I, 1.5 * G),
    Eq(I, 3.0 * J),
    Eq(B - H, -39),
    Eq(D - I, -104),
    Eq(J - F, 59),
    Eq(B - F, -7)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, F, G, H, I, J, K))

# Extract the solution and sort it in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

# Since E is not involved in any equation, we need to assign it one of the remaining numbers
# List of numbers: 3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225
# Find the number that is not used in the solution
used_numbers = set(result)
all_numbers = {3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225}
remaining_number = list(all_numbers - used_numbers)[0]

# Insert the remaining number for E
result.insert(4, remaining_number)

print(result)