from sympy import symbols, Eq, solve

# Define symbols for each letter
A, B, C, E, F, G, H, I = symbols('A B C E F G H I')

# Define the equations based on the constraints
equations = [
    Eq(B - I, -74),
    Eq(C + E, 140),
    Eq(A + G, 39),
    Eq(E + H, 29),
    Eq(B, 2 * A),
    Eq(E + I, 100),
    Eq(A + F, 124),
    Eq(A - H, -6)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, E, F, G, H, I))

# List of possible numbers
numbers = {3, 6, 9, 20, 32, 36, 80, 120, 121}

# Extract the values in alphabetical order
result = [solution[A], solution[B], solution[C], solution[E], solution[F], solution[G], solution[H], solution[I]]

# Determine the value for D
used_numbers = set(result)
D = (numbers - used_numbers).pop()

# Insert D into the correct position in the result list
result.insert(3, D)

# Print the result in the required format
print(f"<<<{result}>>>")