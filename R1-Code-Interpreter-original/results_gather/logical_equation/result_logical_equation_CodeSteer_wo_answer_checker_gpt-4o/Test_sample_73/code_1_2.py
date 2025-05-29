from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(M, 3 * L),
    Eq(F, 2.5 * H),
    Eq(E - L, 9),
    Eq(A + I, 60),
    Eq(E + J, 17),
    Eq(G - I, 21),
    Eq(G - A, -9),
    Eq(E + G, 46),
    Eq(H + I, 17)
]

# Solve the system of equations
solution = solve(equations, (A, E, F, G, H, I, J, L, M))

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Extract the solution values
solution_values = {var: solution[var] for var in solution}

# Assign the numbers to the variables based on the solution
assigned_values = {}
for var, value in solution_values.items():
    for number in numbers:
        if value == number:
            assigned_values[var] = number
            numbers.remove(number)
            break

# Output the values in alphabetical order
result = [assigned_values.get(letter, None) for letter in [A, B, C, D, E, F, G, H, I, J, K, L, M]]
print(f"<<<{result}>>>")