from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if F != 4.0 * C:  # F = 4.0C
        return False
    if J - A != -48:  # J - A = -48
        return False
    if K - B != -75:  # K - B = -75
        return False
    if J <= C:  # J > C
        return False
    if B != 3.0 * A:  # B = 3.0A
        return False
    if B + C != 234:  # B + C = 234
        return False
    if F + J != 63:  # F + J = 63
        return False
    if H + K != 271:  # H + K = 271
        return False
    if A - I != 59:  # A - I = 59
        return False
    if C != 3.0 * E:  # C = 3.0E
        return False
    
    return True

# Try all permutations
for p in permutations(numbers):
    if check_conditions(p):
        print(list(p))
        break