def is_valid_partial(assignment, index):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check constraints that can be evaluated with the current partial assignment
    if index >= 4 and F is not None and E is not None and F != 2 * E:
        return False
    if index >= 1 and B is not None and E is not None and B <= E:
        return False
    if index >= 0 and A is not None and H is not None and A - H != 5:
        return False
    if index >= 3 and D is not None and L is not None and D + L != 48:
        return False
    if index >= 1 and B is not None and J is not None and B <= J:
        return False
    if index >= 1 and B is not None and J is not None and B + J != 16:
        return False
    if index >= 0 and A is not None and B is not None and A + B != 22:
        return False
    if index >= 0 and A is not None and E is not None and A != 1.4 * E:
        return False
    if index >= 6 and G is not None and J is not None and G + J != 25:
        return False
    if index >= 2 and C is not None and E is not None and C + E != 21:
        return False
    
    return True

def solve(assignment, numbers, index=0):
    if index == len(assignment):
        print(f"<<<{list(assignment)}>>>")
        return True
    
    for number in numbers:
        if number not in assignment:
            assignment[index] = number
            if is_valid_partial(assignment, index):
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