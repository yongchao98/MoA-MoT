from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if C + H != 124: return False
    if B + F != 305: return False
    if G != 3.0 * J: return False
    if H - G != -45: return False
    if C + J != 137: return False
    if J - A != -11: return False
    if E - J != 59: return False
    if F != 3.0 * E: return False
    if A + H != 30: return False
    if A - C != -94: return False
    
    return True

# Numbers to assign
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break