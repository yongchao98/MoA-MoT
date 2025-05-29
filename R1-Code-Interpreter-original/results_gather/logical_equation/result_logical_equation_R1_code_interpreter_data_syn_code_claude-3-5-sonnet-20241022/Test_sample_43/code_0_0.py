from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K = values
    
    # Check all constraints
    if J - B != 27: return False
    if K - H != 123: return False
    if C - I != -71: return False
    if I - C != 71: return False
    if K != 2.0 * J: return False
    if C - B != -39: return False
    if G - E != -33: return False
    if A + D != 346: return False
    if C != 3.0 * G: return False
    if E - B != -12: return False
    
    return True

# Available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        print(list(perm))
        break