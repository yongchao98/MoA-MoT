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
result = [solution[A], solution[B], solution[C], solution[D], solution[F], solution[G], solution[H], solution[I], solution[J]]

# Since E and K are not involved in any equation, we need to assign them the remaining numbers
# List of numbers: 3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225
# Find the numbers that are not used in the solution
used_numbers = set(result)
all_numbers = {3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225}
remaining_numbers = list(all_numbers - used_numbers)

# Assign the remaining numbers to E and K
result.insert(4, remaining_numbers[0])  # Assign one to E
result.append(remaining_numbers[1])     # Assign the other to K

print(result)