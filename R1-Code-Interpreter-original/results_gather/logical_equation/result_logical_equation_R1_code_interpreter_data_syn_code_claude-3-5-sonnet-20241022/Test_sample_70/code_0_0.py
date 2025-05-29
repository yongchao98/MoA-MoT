from itertools import permutations

def check_constraints(values, mapping):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if M != 2.4 * L:  # M = 2.4L
        return False
    if F != 3.0 * L:  # F = 3.0L
        return False
    if E + I != 8:  # E + I = 8
        return False
    if K - H != -8:  # K - H = -8
        return False
    if H <= I:  # H > I
        return False
    if L + M != 51:  # L + M = 51
        return False
    if C + J != 9:  # C + J = 9
        return False
    if G != 2.0 * I:  # G = 2.0I
        return False
    if A - D != 49:  # A - D = 49
        return False
    if H != 1.5 * K:  # H = 1.5K
        return False
    
    return True

# Available numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm, numbers):
        print(list(perm))
        break