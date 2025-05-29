from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if B + I != 31: return False
    if I != 1.5 * D: return False
    if K - E != -2: return False
    if B + H != 57: return False
    if C + F != 132: return False
    if K - D != -13: return False
    if M != 2.0 * E: return False
    if A + G != 43: return False
    if M <= B: return False
    if D - H != -34: return False
    if F != 4.0 * I: return False
    if I + L != 69: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        break