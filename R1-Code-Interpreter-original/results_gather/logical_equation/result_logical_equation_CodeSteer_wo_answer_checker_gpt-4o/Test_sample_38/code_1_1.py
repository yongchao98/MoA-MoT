from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the constraints
equations = [
    Eq(C + H, 124),
    Eq(B + F, 305),
    Eq(G, 3 * J),
    Eq(H - G, -45),
    Eq(C + J, 137),
    Eq(J - A, -11),
    Eq(E - J, 59),
    Eq(F, 3 * E),
    Eq(A + H, 30),
    Eq(A - C, -94)
]

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

# Output the result
print(f"<<<{result}>>>")