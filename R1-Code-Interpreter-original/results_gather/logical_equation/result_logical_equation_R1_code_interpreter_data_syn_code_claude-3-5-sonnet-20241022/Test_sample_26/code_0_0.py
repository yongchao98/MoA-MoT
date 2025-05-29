from itertools import permutations

def check_constraints(p):
    A, B, C, D, E, F, G, H, I = p
    
    # Check all constraints
    if B - I != -74: return False
    if C + E != 140: return False
    if A + G != 39: return False
    if E + H != 29: return False
    if B != 2.0 * A: return False
    if E + I != 100: return False
    if A + F != 124: return False
    if A - H != -6: return False
    
    return True

# Numbers that can be assigned
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all possible permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break