from sympy import symbols, Eq, solve, S

# Define symbols for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K', integer=True)

# Define the equations based on the given conditions
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
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K), dict=True, domain=S.Integers)

# Filter solutions to match the possible values
valid_solution = None
for sol in solution:
    if all(sol[var] in possible_values for var in sol):
        valid_solution = sol
        break

# Extract the values in alphabetical order
if valid_solution:
    result = [valid_solution[A], valid_solution[B], valid_solution[C], valid_solution[D], 
              valid_solution[E], valid_solution[F], valid_solution[G], valid_solution[H], 
              valid_solution[I], valid_solution[J], valid_solution[K]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")