from itertools import permutations

def check_equations(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all equations
    if G != 3.0 * C:  # G = 3.0C
        return False
    if I != 4.0 * G:  # I = 4.0G
        return False
    if A - G != 112:  # A - G = 112
        return False
    if K != 3.0 * J:  # K = 3.0J
        return False
    if B + D != 64:   # B + D = 64
        return False
    if E - B != 11:   # E - B = 11
        return False
    if F != 2.0 * J:  # F = 2.0J
        return False
    if F - C != 147:  # F - C = 147
        return False
    if E - I != -9:   # E - I = -9
        return False
    if A + F != 271:  # A + F = 271
        return False
    
    return True

numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_equations(perm):
        print(list(perm))
        break