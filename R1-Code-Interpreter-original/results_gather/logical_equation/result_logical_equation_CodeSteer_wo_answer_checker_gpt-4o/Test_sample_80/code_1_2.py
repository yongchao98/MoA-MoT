from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Constraints as functions
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        D == 3 * L and
        F + L == 17 and
        I + J == 33 and
        I - D == -17 and
        I - E == 21 and
        D - C == 29 and
        A == 4 * B and
        F - C == -14 and
        L == 3 * J and
        L - D == -30 and
        F - B == -22
    )

# Backtracking function
def backtrack(assignment, index):
    if index == len(assignment):
        if satisfies_constraints(assignment):
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

# Initialize assignment list
assignment = [None] * 13

# Find the solution using backtracking
solution = backtrack(assignment, 0)

# Print the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")