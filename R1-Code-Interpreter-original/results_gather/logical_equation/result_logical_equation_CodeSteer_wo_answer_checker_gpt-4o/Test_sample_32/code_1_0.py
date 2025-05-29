from sympy import symbols, Eq, solve

# Define the symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given inequalities
equations = [
    Eq(D - A, 11),
    Eq(H + I, 228),
    Eq(G, 4 * B),
    Eq(I - F, -118),
    Eq(H - B, 216),
    Eq(A + G, 52),
    Eq(F + G, 157),
    Eq(E, 3 * A),
    Eq(G + I, 39),
    Eq(C - H, -75)
]

# Define the set of possible values
possible_values = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Filter solutions to ensure they are within the possible values
valid_solution = None
for sol in solution:
    if all(value in possible_values for value in sol.values()):
        valid_solution = sol
        break

# Extract the solution in alphabetical order if a valid solution is found
if valid_solution:
    result = [valid_solution.get(var, None) for var in (A, B, C, D, E, F, G, H, I, J, K)]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")