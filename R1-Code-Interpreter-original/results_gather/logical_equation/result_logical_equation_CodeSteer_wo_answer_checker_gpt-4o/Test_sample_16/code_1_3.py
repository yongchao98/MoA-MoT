from sympy import symbols, Eq, solve

# Define the variables
A, C, D, E, F, G, H, I = symbols('A C D E F G H I')

# Define the equations based on the constraints
equations = [
    Eq(C, 3.0 * D),
    Eq(C, 1.5 * F),
    Eq(D - C, -6),
    Eq(F, 2.0 * D),
    Eq(D + E, 124),
    Eq(A + H, 112),
    Eq(I, 1.5 * A),
    Eq(D + G, 23)
]

# Solve the equations
solution = solve(equations, (A, C, D, E, F, G, H, I), dict=True)

# Given numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Extract the solution and verify it
if solution:
    # Since solve returns a list of solutions, we need to iterate over them
    for sol in solution:
        # Extract the values in alphabetical order
        result = [sol[A], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I]]
        
        # Check if the solution uses all the given numbers
        if sorted(result) == sorted(numbers):
            print(f"<<<{result}>>>")
            break
    else:
        print("No valid solution found that matches the given numbers.")
else:
    print("No solution found.")