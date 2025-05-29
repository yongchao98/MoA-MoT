from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if C + H != 52: return False
    if G + I != 13: return False
    if abs(H - 2.4*E) > 0.0001: return False  # Using small epsilon for float comparison
    if F + J != 57: return False
    if A - I != 25: return False
    if J - I != 4: return False
    if D + G != 34: return False
    if F + L != 95: return False
    if abs(A - 4.0*J) > 0.0001: return False
    if F + G != 60: return False
    if abs(J - 1.4*K) > 0.0001: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break