from itertools import permutations

# Define the numbers and the letters
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Define a function to check if a partial assignment is valid
def is_valid(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        A - E == 4 and
        I == 4 * J and
        D - B == -14 and
        H + M == 26 and
        L == 3.2 * A and
        L == 1.6 * F and
        F + L == 26 and
        H - I == -4 and
        D == 1.5 * H
    )

# Use backtracking to find a valid assignment
def backtrack(assignment, numbers):
    if len(assignment) == 13:
        if is_valid(assignment):
            return assignment
        return None

    for number in numbers:
        if number not in assignment:
            result = backtrack(assignment + [number], numbers)
            if result is not None:
                return result
    return None

# Start the backtracking process
solution = backtrack([], numbers)

# Print the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")