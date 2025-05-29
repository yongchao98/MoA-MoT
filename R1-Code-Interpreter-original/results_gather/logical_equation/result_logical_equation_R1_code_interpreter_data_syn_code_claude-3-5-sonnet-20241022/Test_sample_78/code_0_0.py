from itertools import permutations

def check_constraints(values, mapping):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if H != 2.4 * D:  # H = 2.4D
        return False
    if G - A != -9:  # G - A = -9
        return False
    if E != 4.0 * G:  # E = 4.0G
        return False
    if F - G != 8:  # F - G = 8
        return False
    if B + E != 78:  # B + E = 78
        return False
    if M - K != 40:  # M - K = 40
        return False
    if D + M != 55:  # D + M = 55
        return False
    if abs(A - 3.2 * K) > 0.0001:  # A = 3.2K (using small epsilon for float comparison)
        return False
    if I - C != 94:  # I - C = 94
        return False
    if G + L != 43:  # G + L = 43
        return False
    if L - A != 20:  # L - A = 20
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