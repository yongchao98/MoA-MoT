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
    
    # Known values from the solution
    known_values = {var: solution[var] for var in solution}
    
    # Determine remaining values
    all_values = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}
    used_values = set(known_values.values())
    remaining_values = list(all_values - used_values)
    
    # Debugging: Print known values and remaining values
    print("Known values:", known_values)
    print("Remaining values:", remaining_values)
    
    # Assign remaining values to unsolved variables
    unsolved_vars = [A, B, J]
    for var, value in zip(unsolved_vars, remaining_values):
        known_values[var] = value
    
    # Output the values in alphabetical order
    result = [known_values[var] for var in [A, B, C, D, E, F, G, H, I, J, K]]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")