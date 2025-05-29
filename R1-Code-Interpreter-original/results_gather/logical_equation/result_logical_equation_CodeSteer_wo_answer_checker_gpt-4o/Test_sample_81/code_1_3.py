from sympy import symbols, Eq
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations based on the given constraints
equations = [
    Eq(F + G, 46),
    Eq(J - M, 80),
    Eq(J - F, 86),
    Eq(H - L, 26),
    Eq(B + M, 31),
    Eq(B, 1.5 * F),
    Eq(C, 4.0 * A),
    Eq(L, 1.5 * M),
    Eq(G, 3.6 * F),
    Eq(D, 1.5 * E),
    Eq(L - I, -21)
]

# Define the possible values for each letter
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Function to check if the current assignment satisfies all equations
def satisfies_constraints(assignment):
    return all(eq.subs(assignment) for eq in equations)

# Backtracking function
def backtrack(assignment, variables, values):
    if len(assignment) == len(variables):
        if satisfies_constraints(assignment):
            return assignment
        return None
    
    var = variables[len(assignment)]
    for value in values:
        if value not in assignment.values():
            assignment[var] = value
            result = backtrack(assignment, variables, values)
            if result is not None:
                return result
            del assignment[var]
    return None

# Solve the problem using backtracking
variables = [A, B, C, D, E, F, G, H, I, J, K, L, M]
assignment = {}
solution = backtrack(assignment, variables, possible_values)

# Extract the values in alphabetical order and output the result
if solution:
    result = [solution[A], solution[B], solution[C], solution[D],
              solution[E], solution[F], solution[G], solution[H],
              solution[I], solution[J], solution[K], solution[L],
              solution[M]]
    print(f"<<<{result}>>>")
else:
    print("No solution found.")