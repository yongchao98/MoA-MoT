from sympy import symbols, Eq, solve
from itertools import permutations

# Define symbolic variables for each letter
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the equations and inequalities
constraints = [
    Eq(A - G, -13),
    Eq(B + C, 15),
    Eq(H, 4.8 * B),
    Eq(I - B, -3),
    Eq(L - J, -5),
    B > A,
    Eq(I, 2.0 * F),
    Eq(C - H, -14),
    Eq(B - I, 3),
    Eq(A - F, 2)
]

# List of numbers to assign
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Function to check if a partial assignment satisfies all constraints
def satisfies_constraints(assignment):
    # Create a dictionary from the assignment
    assignment_dict = {var: val for var, val in zip([A, B, C, D, E, F, G, H, I, J, K, L, M], assignment)}
    
    # Check each constraint
    for constraint in constraints:
        if not constraint.subs(assignment_dict):
            return False
    return True

# Backtracking function to find a valid assignment
def backtrack(assignment):
    if len(assignment) == 13:
        # If all variables are assigned, check if the assignment satisfies all constraints
        if satisfies_constraints(assignment):
            return assignment
        else:
            return None
    
    # Try assigning each number to the next variable
    for number in numbers:
        if number not in assignment:
            result = backtrack(assignment + [number])
            if result is not None:
                return result
    
    return None

# Start the backtracking process with an empty assignment
solution = backtrack([])

# If a solution is found, print it
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")