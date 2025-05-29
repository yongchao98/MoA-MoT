from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if I - E != -4: return False
    if abs(F - 2.5*K) > 0.01: return False  # Using small epsilon for float comparison
    if abs(I - 1.5*B) > 0.01: return False
    if L - A != 33: return False
    if abs(E - 2.8*H) > 0.01: return False
    if C + M != 111: return False
    if L - C != -60: return False
    if abs(H - 2.0*F) > 0.01: return False
    if F + L != 41: return False
    if K - E != -26: return False
    if A + F != 8: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        break