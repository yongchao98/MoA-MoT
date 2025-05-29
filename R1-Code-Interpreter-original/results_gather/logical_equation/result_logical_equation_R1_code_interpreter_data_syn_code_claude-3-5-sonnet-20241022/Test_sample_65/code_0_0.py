from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if F != 2.0 * E:  # F = 2.0E
        return False
    if B <= E:  # B > E
        return False
    if A - H != 5:  # A - H = 5
        return False
    if D + L != 48:  # D + L = 48
        return False
    if B <= J:  # B > J
        return False
    if B + J != 16:  # B + J = 16
        return False
    if A + B != 22:  # A + B = 22
        return False
    if abs(A - 1.4 * E) > 0.0001:  # A = 1.4E
        return False
    if G + J != 25:  # G + J = 25
        return False
    if C + E != 21:  # C + E = 21
        return False
    
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        print(list(perm))
        break