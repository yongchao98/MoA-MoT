from itertools import permutations

def check_constraints(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if C - E != -21: return False
    if J != 4.0 * C: return False
    if H + K != 23: return False
    if D + L != 17: return False
    if C - A != -26: return False
    if G + M != 8: return False
    if F + G != 13: return False
    if H - G != 4: return False
    if A - C != 26: return False
    if I - J != -68: return False
    if F + M != 15: return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break