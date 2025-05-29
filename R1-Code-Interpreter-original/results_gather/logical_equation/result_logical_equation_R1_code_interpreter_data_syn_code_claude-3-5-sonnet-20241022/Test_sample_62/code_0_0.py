from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if G != 3.6 * J:  # G = 3.6J
        return False
    if C != 3.0 * E:  # C = 3.0E
        return False
    if E + I != 18:  # E + I = 18
        return False
    if K - C != 5:  # K - C = 5
        return False
    if L - B != 17:  # L - B = 17
        return False
    if B + H != 12:  # B + H = 12
        return False
    if B + C != 52:  # B + C = 52
        return False
    if D - E != -13:  # D - E = -13
        return False
    if B - A != -9:  # B - A = -9
        return False
    if A - J != 6:  # A - J = 6
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break