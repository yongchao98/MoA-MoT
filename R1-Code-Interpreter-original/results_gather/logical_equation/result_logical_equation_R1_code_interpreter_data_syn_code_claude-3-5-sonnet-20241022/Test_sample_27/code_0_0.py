from itertools import permutations

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

def check_constraints(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if not abs(C - 3.0 * H) < 0.001:  # C = 3.0H
        return False
    if A - I != 16:  # A - I = 16
        return False
    if I - H != 17:  # I - H = 17
        return False
    if not abs(B - 4.0 * I) < 0.001:  # B = 4.0I
        return False
    if E + F != 127:  # E + F = 127
        return False
    if F - G != -114:  # F - G = -114
        return False
    if not abs(A - 4.0 * C) < 0.001:  # A = 4.0C
        return False
    if not abs(B - 2.5 * D) < 0.001:  # B = 2.5D
        return False
    
    return True

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break