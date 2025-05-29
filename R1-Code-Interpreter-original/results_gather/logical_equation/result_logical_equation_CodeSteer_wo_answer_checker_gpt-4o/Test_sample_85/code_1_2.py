def is_valid(assignment):
    # Unpack the assignment
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    return (
        E > C and
        G == 4.8 * M and
        L - C == 12 and
        D + H == 66 and
        L - D == -35 and
        E == 3.6 * B and
        B - D == -40 and
        H - G == -8 and
        C == 1.5 * K and
        H == 1.6 * B and
        J - E == 9
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
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Initialize assignment with None
assignment = [None] * 13

# Find the solution using backtracking
solution = backtrack(assignment, numbers, 0)

# Print the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")