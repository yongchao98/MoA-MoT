from sympy import symbols, Eq, solve
from itertools import permutations

# Define the variables
A, B, C, D, E, F, G, H, I, J, K, L, M = symbols('A B C D E F G H I J K L M')

# Define the constraints
constraints = [
    Eq(D, 3 * I),
    Eq(C + H, 20),
    A > F,
    Eq(F + G, 38),
    Eq(L - G, -26),
    Eq(C + F, 17),
    Eq(E + L, 60),
    Eq(M, 4.5 * L),
    Eq(A - K, -17),
    Eq(D + H, 8)
]

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Function to check if a partial assignment satisfies the constraints
def check_constraints(assignment):
    subs = {var: val for var, val in zip([A, B, C, D, E, F, G, H, I, J, K, L, M], assignment)}
    return all(constraint.subs(subs) for constraint in constraints)

# Backtracking function with constraint propagation
def backtrack(assignment, index):
    if index == len(numbers):
        if check_constraints(assignment):
            return assignment
        return None
    
    for num in numbers:
        if num not in assignment:
            assignment[index] = num
            result = backtrack(assignment, index + 1)
            if result is not None:
                return result
            assignment[index] = None
    
    return None

# Initialize assignment
assignment = [None] * 13

# Find solution using backtracking
solution = backtrack(assignment, 0)

# Output the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")