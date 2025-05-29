def is_valid_partial(assignment, index):
    # Unpack the assignment
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check constraints that can be evaluated with the current partial assignment
    if index >= 4 and E <= C:
        return False
    if index >= 6 and G != 4.8 * M:
        return False
    if index >= 11 and L - C != 12:
        return False
    if index >= 7 and D + H != 66:
        return False
    if index >= 11 and L - D != -35:
        return False
    if index >= 4 and E != 3.6 * B:
        return False
    if index >= 3 and B - D != -40:
        return False
    if index >= 7 and H - G != -8:
        return False
    if index >= 2 and C != 1.5 * K:
        return False
    if index >= 7 and H != 1.6 * B:
        return False
    if index >= 9 and J - E != 9:
        return False
    
    return True

def backtrack(assignment, numbers, index):
    if index == len(assignment):
        return assignment if is_valid_partial(assignment, index) else None
    
    for number in numbers:
        if number not in assignment:
            assignment[index] = number
            if is_valid_partial(assignment, index):
                result = backtrack(assignment, numbers, index + 1)
                if result is not None:
                    return result
            assignment[index] = None
    
    return None

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Initialize assignment with None
assignment = [None] * 13

# Find the solution using backtracking with forward checking
solution = backtrack(assignment, numbers, 0)

# Print the solution
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")