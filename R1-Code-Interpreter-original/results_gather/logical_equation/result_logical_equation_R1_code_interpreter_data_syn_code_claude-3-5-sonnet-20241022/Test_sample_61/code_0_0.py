from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if C + M != 43: return False
    if abs(C - 1.5*L) > 0.001: return False  # Using small epsilon for float comparison
    if abs(L - 2.0*J) > 0.001: return False
    if A + K != 81: return False
    if C - B != -35: return False
    if H - K != -20: return False
    if abs(F - 4.8*J) > 0.001: return False
    if K - M != 8: return False
    if D + E != 4: return False
    if abs(K - 2.4*C) > 0.001: return False
    
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break