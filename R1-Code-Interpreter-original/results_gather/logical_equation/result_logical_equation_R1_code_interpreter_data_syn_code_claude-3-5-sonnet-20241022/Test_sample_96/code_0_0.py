from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if L != 1.4 * K:  # L = 1.4K
        return False
    if B + C != 73:  # B + C = 73
        return False
    if J - G != -86:  # J - G = -86
        return False
    if D + H != 53:  # D + H = 53
        return False
    if B + K != 50:  # B + K = 50
        return False
    if L - A != -29:  # L - A = -29
        return False
    if F - I != -13:  # F - I = -13
        return False
    if J != 2.0 * K:  # J = 2.0K
        return False
    if C + M != 52:  # C + M = 52
        return False
    if I <= L:  # I > L
        return False
    if D + M != 74:  # D + M = 74
        return False
    if M - B != -21:  # M - B = -21
        return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        print(list(perm))
        break