from sympy import symbols, Eq, nonlinsolve

# Define the variables
A, B, C, D, E, F, G, H, I = symbols('A B C D E F G H I')

# Define the equations based on the given conditions
equations = [
    Eq(G, 1.5 * H),
    Eq(E - D, -41),
    Eq(F, 1.5 * E),
    Eq(D + G, 130),
    Eq(F + H, 126),
    Eq(F + I, 123),
    Eq(E, 2.5 * B),
    Eq(D - E, 41)
]

# Solve the system of equations
solution_set = nonlinsolve(equations, (A, B, C, D, E, F, G, H, I))

# Extract the solution
for solution in solution_set:
    # Convert the solution to a list and print it in the required format
    result = list(solution)
    print(f"<<<{result}>>>")
    break  # Assuming we only need the first valid solution