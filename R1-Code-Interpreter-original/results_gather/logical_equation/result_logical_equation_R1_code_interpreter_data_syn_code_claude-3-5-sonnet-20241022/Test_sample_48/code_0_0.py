from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K = values
    
    # Check all constraints
    if G + K != 23: return False  # 1
    if I != 3.0 * C: return False  # 2
    if D - F != 7: return False  # 3
    if D - B != -23: return False  # 4
    if J != 2.4 * K: return False  # 5
    if J <= B: return False  # 6
    if D + I != 241: return False  # 7
    if E - G != 118: return False  # 8
    if B - K != 19: return False  # 9
    if F != 3.0 * G: return False  # 10
    if C - J != 27: return False  # 11
    
    return True

# Available numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        print(list(perm))
        break