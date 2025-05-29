from itertools import permutations

def check_constraints(p):
    A, B, C, D, E, F, G, H, I = p
    
    # Check all constraints
    if E + F != 83: return False
    if A - H != -112: return False
    if G != 4.0 * C: return False
    if B + D != 63: return False
    if D != 3.0 * A: return False
    if G + H != 125: return False
    if B - E != 33: return False
    
    return True

# Numbers that can be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break