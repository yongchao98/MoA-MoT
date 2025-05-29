from sympy import symbols, Eq, solve

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the constraints
equations = [
    Eq(E, 1.5 * K),
    Eq(D + I, 155),
    Eq(K, 2.0 * D),
    Eq(E, 3.0 * D),
    Eq(D - C, -46),
    Eq(F + H, 51),
    Eq(H + I, 83),
    Eq(F, 3.0 * G),
    Eq(E + K, 375),
    Eq(F - D, -27)
]

# Solve the system of equations
solution = solve(equations, (C, D, E, F, G, H, I, K), dict=True)

# Check if a solution was found
if solution:
    # Extract the first solution (if multiple solutions exist)
    solution = solution[0]
    
    # Assign known values
    C = solution[C]
    D = solution[D]
    E = solution[E]
    F = solution[F]
    G = solution[G]
    H = solution[H]
    I = solution[I]
    K = solution[K]
    
    # Determine remaining values
    all_values = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}
    used_values = set(solution.values())
    remaining_values = list(all_values - used_values)
    
    # Debugging: Print remaining values
    print("Remaining values:", remaining_values)
    
    # Assign remaining values to A, B, J
    # Ensure the number of remaining values matches the number of variables
    if len(remaining_values) == 3:
        A, B, J = remaining_values
    else:
        raise ValueError("Unexpected number of remaining values.")
    
    # Output the values in alphabetical order
    result = [A, B, C, D, E, F, G, H, I, J, K]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")