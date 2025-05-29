from sympy import symbols, Eq, solve

# Define the symbols for the relevant variables
A, B, C, D, E, F, H, I = symbols('A B C D E F H I')

# Define the equations based on the given conditions
equations = [
    Eq(A + C, 127),
    Eq(D - E, -77),
    Eq(A - B, 89),
    Eq(B, 1.6 * F),
    Eq(E, 4.0 * F),
    Eq(I, 4.0 * H),
    Eq(H, 1.5 * C)
]

# Define the numbers to be used
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Solve the equations
solution = solve(equations, (A, B, C, D, E, F, H, I))

# Extract the solution and check against the numbers
result = [solution[A], solution[B], solution[C], solution[D], solution[E], solution[F], solution[H], solution[I]]

# Check if the solution matches the available numbers
if sorted(result) == sorted(numbers):
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")