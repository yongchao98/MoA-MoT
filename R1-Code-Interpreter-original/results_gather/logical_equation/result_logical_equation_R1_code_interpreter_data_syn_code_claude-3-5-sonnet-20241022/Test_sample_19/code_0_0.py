from itertools import permutations

def check_constraints(values, numbers):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all constraints
    if A + F != 141: return False
    if A + G != 124: return False
    if I - C != 26: return False
    if H - G != 33: return False
    if C != 2.0 * G: return False
    if D - H != -27: return False
    if D - G != 6: return False
    if H - E != -84: return False
    
    return True

# Numbers that can be assigned
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm, numbers):
        print(list(perm))
        break