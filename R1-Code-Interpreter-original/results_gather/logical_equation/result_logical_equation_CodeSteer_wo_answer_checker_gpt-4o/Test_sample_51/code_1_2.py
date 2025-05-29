from sympy import symbols, Eq, nonlinsolve, S

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K', integer=True)

# Define the equations and inequalities based on the given conditions
equations = [
    Eq(K, 2 * E),
    Eq(E - C, 72),
    Eq(J - G, -205),
    Eq(A + F, 128),
    Eq(C + F, 51),
    Eq(G, 1.5 * K),
    Eq(B + C, 19),
    Eq(H, 3 * C),
    Eq(F, 2.4 * J),
    Eq(G + J, 245),
    K > H
]

# Define the possible values for each letter
possible_values = {3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225}

# Solve the system of equations
solution_set = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Filter solutions to match the possible values
valid_solution = None
for sol in solution_set:
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