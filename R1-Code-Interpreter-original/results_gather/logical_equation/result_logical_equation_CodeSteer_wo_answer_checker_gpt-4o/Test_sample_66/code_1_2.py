from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Constraints as functions
def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    return (
        D == 3.0 * I and
        C + H == 20 and
        A > F and
        F + G == 38 and
        L - G == -26 and
        C + F == 17 and
        E + L == 60 and
        M == 4.5 * L and
        A - K == -17 and
        D + H == 8
    )

# Backtracking function
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