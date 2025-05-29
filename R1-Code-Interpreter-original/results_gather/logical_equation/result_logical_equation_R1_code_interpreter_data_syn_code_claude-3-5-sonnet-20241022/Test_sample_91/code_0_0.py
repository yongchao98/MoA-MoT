from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if not abs(I - 2.4*F) < 0.1: return False  # I = 2.4F
    if not abs(D - 4.8*C) < 0.1: return False  # D = 4.8C
    if not abs(D - 1.6*F) < 0.1: return False  # D = 1.6F
    if not abs(B - 1.6*L) < 0.1: return False  # B = 1.6L
    if not abs(M - 3.5*E) < 0.1: return False  # M = 3.5E
    if not G - J == -17: return False          # G - J = -17
    if not H - B == 34: return False           # H - B = 34
    if not abs(B - 3.2*C) < 0.1: return False  # B = 3.2C
    if not H > E: return False                 # H > E
    if not H + L == 60: return False          # H + L = 60
    if not I - C == 31: return False          # I - C = 31
    if not H > A: return False                # H > A
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break