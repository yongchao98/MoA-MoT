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

# Solve the system of equations
solution = solve(equations, (A, B, C, D, E, F, G, H, I, J, K))

# Extract the solution in alphabetical order
result = [solution[A], solution[B], solution[C], solution[D], solution[E], 
          solution[F], solution[G], solution[H], solution[I], solution[J], solution[K]]

# Print the result in the required format
print(f"<<<{result}>>>")