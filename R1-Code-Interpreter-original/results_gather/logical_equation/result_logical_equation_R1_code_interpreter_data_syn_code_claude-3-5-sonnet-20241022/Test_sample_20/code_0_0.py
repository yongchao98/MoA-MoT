from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I = [values[letters.index(x)] for x in 'ABCDEFGHI']
    
    # Check all constraints
    if E + G != 127: return False
    if H != 4.0 * D: return False
    if A + D != 140: return False
    if B != 1.6 * D: return False
    if A - G != -1: return False
    if G <= B: return False
    if C + D != 29: return False
    if G - C != 112: return False
    
    return True

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = list('ABCDEFGHI')

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        # Create result in alphabetical order
        result = list(perm)
        print(result)
        break