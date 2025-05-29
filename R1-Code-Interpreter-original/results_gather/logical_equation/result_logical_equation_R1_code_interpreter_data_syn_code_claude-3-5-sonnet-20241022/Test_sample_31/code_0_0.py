from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if C - I != 13: return False
    if G - I != 147: return False
    if abs(A - 3.0 * D) > 0.1: return False  # Using small epsilon for float comparison
    if G + J != 186: return False
    if D + H != 123: return False
    if F - A != -216: return False
    if H + J != 84: return False
    if F + G != 159: return False
    if A + B != 346: return False
    if abs(F - 3.0 * I) > 0.1: return False  # Using small epsilon for float comparison
    
    return True

# Available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break