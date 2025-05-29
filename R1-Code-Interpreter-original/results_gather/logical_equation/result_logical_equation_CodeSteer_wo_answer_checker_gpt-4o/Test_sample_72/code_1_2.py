def is_valid(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        F + J == 86 and
        B == 4.5 * A and
        A + K == 15 and
        M == 3.0 * L and
        H == 1.5 * A and
        J > D and
        A - M == 7 and
        B + E == 52 and
        A > G and
        H > L
    )

def backtrack(assignment, numbers, index):
    if index == len(assignment):
        if is_valid(assignment):
            return assignment
        return None

    for number in numbers:
        if number not in assignment:
            assignment[index] = number
            result = backtrack(assignment, numbers, index + 1)
            if result is not None:
                return result
            assignment[index] = None

    return None

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Initialize assignment with None
assignment = [None] * 13

# Solve the CSP
solution = backtrack(assignment, numbers, 0)

# Output the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")