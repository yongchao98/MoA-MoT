from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if A + E != 107: return False
    if D + H != 137: return False
    if H <= C: return False
    if A + B != 63: return False
    if E - F != 79: return False
    if B + E != 116: return False
    if B - D != -85: return False
    
    return True

# Numbers to assign
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break