from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if K - G != 25: return False
    if F != 3.5 * M: return False
    if B + E != 141: return False
    if J + K != 38: return False
    if B + I != 50: return False
    if A - G != 33: return False
    if I + M != 7: return False
    if A - H != -14: return False
    if M - H != -48: return False
    if C - E != -72: return False
    if F != 1.4 * I: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break