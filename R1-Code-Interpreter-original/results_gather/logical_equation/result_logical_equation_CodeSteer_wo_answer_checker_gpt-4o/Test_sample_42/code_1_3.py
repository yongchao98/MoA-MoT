from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given inequalities
equations = [
    Eq(H + J, 130),
    Eq(C + I, 177),
    Eq(C + G, 75),
    Eq(G, 3 * F),
    Eq(B - E, 150),
    Eq(J - C, -18),
    Eq(E + G, 123),
    Eq(D + F, 96),
    Eq(E + K, 78),
    Eq(H + K, 124)
]

# Solve the system of equations
solutions = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Possible values for the letters
possible_values = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}

# Iterate over the solutions
for solution in solutions:
    # Extract the values using symbolic variables
    values = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]
    
    # Check if the values match the possible values
    if set(values) == possible_values:
        print(f"<<<{values}>>>")
        break