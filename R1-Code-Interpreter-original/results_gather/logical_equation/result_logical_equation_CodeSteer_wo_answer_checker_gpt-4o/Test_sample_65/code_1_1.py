def is_valid(assignment):
    # Unpack the assignment
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    return (
        F == 2 * E and
        B > E and
        A - H == 5 and
        D + L == 48 and
        B > J and
        B + J == 16 and
        A + B == 22 and
        A == 1.4 * E and
        G + J == 25 and
        C + E == 21
    )

def solve(assignment, numbers, index=0):
    if index == len(assignment):
        if is_valid(assignment):
            print(f"<<<{list(assignment)}>>>")
            return True
        return False
    
    for number in numbers:
        if number not in assignment:
            assignment[index] = number
            if solve(assignment, numbers, index + 1):
                return True
            assignment[index] = None
    
    return False

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Initialize assignment with None
assignment = [None] * 13

# Solve the problem
solve(assignment, numbers)