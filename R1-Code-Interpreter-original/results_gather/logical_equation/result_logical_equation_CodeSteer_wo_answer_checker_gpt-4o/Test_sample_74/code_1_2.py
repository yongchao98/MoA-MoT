from sympy import symbols, Eq, nonlinsolve, S

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M', integer=True)

# Define the equations
equations = [
    Eq(E + F, 8),
    Eq(C - D, 14),
    Eq(A + M, 35),
    Eq(D + K, 3),
    Eq(G + J, 60),
    Eq(H + I, 60),
    Eq(A - K, 5),
    Eq(G - H, 26),
    Eq(G - F, 45),
    B > F
]

# Define the domain of possible values
possible_values = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}

# Use nonlinsolve to find solutions
solutions = nonlinsolve(equations, [A, B, C, D, E, F, G, H, I, J, K, L, M])

# Filter solutions to ensure they meet the constraints
valid_solution = None
for sol in solutions:
    if all(val in possible_values for val in sol):
        valid_solution = sol
        break

# Extract the solution in alphabetical order if a valid solution is found
if valid_solution:
    result = list(valid_solution)
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")