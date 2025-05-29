from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if D - A != 11: return False
    if H + I != 228: return False
    if G != 4.0 * B: return False
    if I - F != -118: return False
    if H - B != 216: return False
    if A + G != 52: return False
    if F + G != 157: return False
    if E != 3.0 * A: return False
    if G + I != 39: return False
    if C - H != -75: return False
    
    return True

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break