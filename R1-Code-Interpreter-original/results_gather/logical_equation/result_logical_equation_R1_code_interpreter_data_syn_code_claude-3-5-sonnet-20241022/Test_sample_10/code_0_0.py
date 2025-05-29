from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if C - D != -94:  # Condition 1
        return False
    if C - F != 18:   # Condition 2
        return False
    if B <= E:        # Condition 3
        return False
    if A <= E:        # Condition 4
        return False
    if G != 4 * E:    # Condition 5
        return False
    if B <= F:        # Condition 6
        return False
    if C != 3 * F:    # Condition 7
        return False
    
    return True

# Available numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break