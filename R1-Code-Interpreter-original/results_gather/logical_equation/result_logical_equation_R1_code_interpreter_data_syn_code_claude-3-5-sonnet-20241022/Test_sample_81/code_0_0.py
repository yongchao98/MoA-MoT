from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if F + G != 46: return False
    if J - M != 80: return False
    if J - F != 86: return False
    if H - L != 26: return False
    if B + M != 31: return False
    if abs(B - 1.5 * F) > 0.0001: return False
    if abs(C - 4.0 * A) > 0.0001: return False
    if abs(L - 1.5 * M) > 0.0001: return False
    if abs(G - 3.6 * F) > 0.0001: return False
    if abs(D - 1.5 * E) > 0.0001: return False
    if L - I != -21: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break