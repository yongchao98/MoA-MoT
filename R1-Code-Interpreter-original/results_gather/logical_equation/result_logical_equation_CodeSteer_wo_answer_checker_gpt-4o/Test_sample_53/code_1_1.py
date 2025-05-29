from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
equations = [
    Eq(F - K, 36),
    Eq(B - K, 111),
    Eq(B + E, 230),
    Eq(H + J, 25),
    Eq(E - A, 77),
    Eq(H + I, 36),
    Eq(C + E, 305),
    Eq(C - I, 205),
    Eq(C + J, 234),
    Eq(G, 3 * H),
    Eq(A + K, 42)
]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# The set of numbers
numbers = {3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225}

# Filter the solution to ensure all values are in the given set
valid_solution = {var: val for var, val in solution.items() if val in numbers}

# Check if we have a valid solution for all variables
if len(valid_solution) == 11:
    # Extract the values in alphabetical order
    result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D],
              valid_solution[E], valid_solution[F], valid_solution[G], valid_solution[H],
              valid_solution[I], valid_solution[J], valid_solution[K]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")