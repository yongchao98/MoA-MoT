from itertools import permutations

def check_constraints(p):
    # Map numbers to letters for easier reading
    A, B, C, D, E, F, G, H, I, J, K, L, M = p
    
    # Check all constraints
    if G + I != 55: return False        # 1
    if E <= D: return False             # 2
    if G + L != 57: return False        # 3
    if G - E != 5: return False         # 4
    if C + H != 27: return False        # 5
    if A - I != 5: return False         # 6
    if F - I != 91: return False        # 7
    if B - C != -9: return False        # 8
    if K <= C: return False             # 9
    if H + L != 10: return False        # 10
    if J - D != 14: return False        # 11
    if B - F != -81: return False       # 12
    
    return True

# Numbers that can be assigned
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break