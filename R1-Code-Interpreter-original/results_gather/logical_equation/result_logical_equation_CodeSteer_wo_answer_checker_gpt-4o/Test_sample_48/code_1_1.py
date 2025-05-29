from sympy import symbols, Eq, solve

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations based on the given constraints
equations = [
    Eq(G + K, 23),
    Eq(I, 3 * C),
    Eq(D - F, 7),
    Eq(D - B, -23),
    Eq(J, 2.4 * K),
    Eq(D + I, 241),
    Eq(E - G, 118),
    Eq(B - K, 19),
    Eq(F, 3 * G),
    Eq(C - J, 27)
]

# Define the inequality
inequality = J > B

# Solve the system of equations
solution = solve(equations + [inequality], (A, B, C, D, E, F, G, H, I, J, K), dict=True)

# Since we have a specific set of numbers, we need to map the solution to these numbers
# The numbers are: 3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Extract the solution and map it to the numbers
if solution:
    # Extract the first solution (if multiple solutions exist)
    sol = solution[0]
    # Assign the numbers to the variables based on the solution
    result = [sol[A], sol[B], sol[C], sol[D], sol[E], sol[F], sol[G], sol[H], sol[I], sol[J], sol[K]]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")