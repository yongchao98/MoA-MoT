from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if B - F != -21: return False
    if F + L != 73: return False
    if F - A != -22: return False
    if K - L != -42: return False
    if B + C != 43: return False
    if H + L != 60: return False
    if H + M != 17: return False
    if abs(C - 2.4*H) > 0.01: return False  # Using small epsilon for float comparison
    if D + K != 4: return False
    if abs(G - 2.5*M) > 0.01: return False  # Using small epsilon for float comparison
    
    return True

# Available numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break