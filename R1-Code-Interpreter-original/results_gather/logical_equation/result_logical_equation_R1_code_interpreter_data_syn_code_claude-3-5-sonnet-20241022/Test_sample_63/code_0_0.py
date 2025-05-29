from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if A - E != 4:
        return False
    if E - A != -4:
        return False
    if abs(I - 4.0*J) > 0.01:  # Using small epsilon for float comparison
        return False
    if D - B != -14:
        return False
    if H + M != 26:
        return False
    if abs(L - 3.2*A) > 0.01:
        return False
    if abs(L - 1.6*F) > 0.01:
        return False
    if F + L != 26:
        return False
    if H - I != -4:
        return False
    if abs(D - 1.5*H) > 0.01:
        return False
    
    return True

# Available numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break