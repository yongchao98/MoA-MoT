from sympy import symbols, Eq, nonlinsolve, S

# Define symbolic variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations and inequalities
equations = [
    Eq(F - C, -141),
    C > D,
    Eq(E + K, 305),
    Eq(K, 3 * A),
    Eq(H, 3 * F),
    Eq(I + J, 157),
    Eq(K, 1.5 * C),
    Eq(K - A, 150),
    Eq(F + G, 25),
    Eq(B + C, 153)
]

# Define the possible values for each letter
possible_values = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Solve the system of equations and inequalities
solutions = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Filter solutions to match the possible values
valid_solution = None
for sol in solutions:
    if all(val in possible_values for val in sol):
        valid_solution = sol
        break

# Extract the values in alphabetical order
if valid_solution:
    result = [valid_solution[0], valid_solution[1], valid_solution[2], valid_solution[3],
              valid_solution[4], valid_solution[5], valid_solution[6], valid_solution[7],
              valid_solution[8], valid_solution[9], valid_solution[10]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")