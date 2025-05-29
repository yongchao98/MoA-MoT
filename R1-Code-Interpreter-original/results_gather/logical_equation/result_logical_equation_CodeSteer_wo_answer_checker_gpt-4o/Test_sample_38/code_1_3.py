from sympy import symbols, Eq, solve

# Define symbolic variables for each letter except D and I, which are unconstrained
A, B, C, E, F, G, H, J = symbols('A B C E F G H J')

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
solution = solve(equations, (A, B, C, E, F, G, H, J))

# List of possible numbers
all_numbers = {3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225}

# Extract used numbers from the solution
used_numbers = {solution[A], solution[B], solution[C], solution[E], solution[F], solution[G], solution[H], solution[J]}

# Determine the remaining numbers for D and I
remaining_numbers = list(all_numbers - used_numbers)

# Assign the remaining numbers to D and I
# Since D and I are unconstrained, we can assign them the remaining numbers
D, I = remaining_numbers

# Extract the solution in alphabetical order, including D and I
result = [solution[A], solution[B], solution[C], D, solution[E], solution[F], solution[G], solution[H], I, solution[J], 225]

# Output the result
print(f"<<<{result}>>>")