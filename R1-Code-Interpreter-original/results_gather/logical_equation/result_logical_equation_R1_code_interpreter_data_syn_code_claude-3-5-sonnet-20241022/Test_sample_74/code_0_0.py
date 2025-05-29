from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if E + F != 8: return False
    if C - D != 14: return False
    if A + M != 35: return False
    if D + K != 3: return False
    if G + J != 60: return False
    if B <= F: return False
    if H + I != 60: return False
    if A - K != 5: return False
    if G - H != 26: return False
    if G - F != 45: return False
    
    return True

# Numbers that can be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break