from itertools import permutations

def check_constraints(values, mapping):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if D != 3.0 * L:  # D = 3.0L
        return False
    if F + L != 17:  # F + L = 17
        return False
    if I + J != 33:  # I + J = 33
        return False
    if I - D != -17:  # I - D = -17
        return False
    if I - E != 21:  # I - E = 21
        return False
    if D - C != 29:  # D - C = 29
        return False
    if A != 4.0 * B:  # A = 4.0B
        return False
    if F - C != -14:  # F - C = -14
        return False
    if L != 3.0 * J:  # L = 3.0J
        return False
    if L - D != -30:  # L - D = -30
        return False
    if F - B != -22:  # F - B = -22
        return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm, numbers):
        print(list(perm))
        break