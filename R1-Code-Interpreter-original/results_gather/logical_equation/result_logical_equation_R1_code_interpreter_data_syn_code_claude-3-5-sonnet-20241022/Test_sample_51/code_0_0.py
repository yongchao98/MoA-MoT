from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K = values
    
    # Check all conditions
    if K != 2.0 * E:  # K = 2.0E
        return False
    if E - C != 72:  # E - C = 72
        return False
    if J - G != -205:  # J - G = -205
        return False
    if A + F != 128:  # A + F = 128
        return False
    if C + F != 51:  # C + F = 51
        return False
    if G != 1.5 * K:  # G = 1.5K
        return False
    if B + C != 19:  # B + C = 19
        return False
    if H != 3.0 * C:  # H = 3.0C
        return False
    if F != 2.4 * J:  # F = 2.4J
        return False
    if G + J != 245:  # G + J = 245
        return False
    if K <= H:  # K > H
        return False
    return True

# Numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break