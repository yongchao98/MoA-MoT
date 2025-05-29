from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(B, 1.4 * G),
    Eq(A, 3.2 * G),
    Eq(J, 3.0 * G),
    Eq(A, 1.6 * K),
    Eq(A + J, 31),
    Eq(A + B, 23),
    Eq(F, 2.4 * K),
    Eq(M - B, 89),
    Eq(I - B, 43),
    Eq(F, 4.8 * G),
    Eq(B, 3.5 * L)
]

# Solve the system of equations
solution = solve(equations, (A, B, F, G, I, J, K, L, M))

# Given numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Extract the solved values
solved_values = {var: solution[var] for var in solution}

# Assign the remaining numbers to C, D, E, H
remaining_vars = [C, D, E, H]
remaining_numbers = [num for num in numbers if num not in solved_values.values()]

# Assign the remaining numbers to the remaining variables
for var, num in zip(remaining_vars, remaining_numbers):
    solved_values[var] = num

# Create a list of values in alphabetical order
values = [solved_values[var] for var in [A, B, C, D, E, F, G, H, I, J, K, L, M]]

# Print the list of values
print(f"<<<{values}>>>")