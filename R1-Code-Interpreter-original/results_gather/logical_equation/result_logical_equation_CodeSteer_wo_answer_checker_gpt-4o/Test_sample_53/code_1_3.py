from sympy import symbols, Eq
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K = symbols('A B C D E F G H I J K')

# Define the equations
equations = [
    Eq(F - K, 36),
    Eq(B - K, 111),
    Eq(B + E, 230),
    Eq(H + J, 25),
    Eq(E - A, 77),
    Eq(H + I, 36),
    Eq(C + E, 305),
    Eq(C - I, 205),
    Eq(C + J, 234),
    Eq(G, 3 * H),
    Eq(A + K, 42)
]

# The set of numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Function to check if current assignment satisfies all equations
def satisfies_all(equations, assignment):
    return all(eq.subs(assignment) for eq in equations)

# Backtracking function
def backtrack(assignment, variables, numbers):
    if len(assignment) == len(variables):
        if satisfies_all(equations, assignment):
            return assignment
        return None

    # Get the next variable to assign
    var = variables[len(assignment)]

    for num in numbers:
        if num not in assignment.values():
            # Try assigning the number to the variable
            assignment[var] = num
            result = backtrack(assignment, variables, numbers)
            if result is not None:
                return result
            # Backtrack
            del assignment[var]

    return None

# Solve the problem using backtracking
variables = [A, B, C, D, E, F, G, H, I, J, K]
assignment = {}
solution = backtrack(assignment, variables, numbers)

# If a solution is found, print it
if solution:
    result = [solution[A], solution[B], solution[C], solution[D],
              solution[E], solution[F], solution[G], solution[H],
              solution[I], solution[J], solution[K]]
    print(f"<<<{result}>>>")
else:
    print("No valid solution found.")