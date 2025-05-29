from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I = [values[letters.index(c)] for c in 'ABCDEFGHI']
    
    # Check all constraints
    if B + I != 156: return False  # 1
    if B + D != 123: return False  # 2
    if D - A != -29: return False  # 3
    if G != 2.0 * D: return False  # 4
    if B != 1.5 * C: return False  # 5
    if A <= D: return False        # 6
    if A - I != -4: return False   # 7
    if C - F != 71: return False   # 8
    
    return True

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = list('ABCDEFGHI')

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        # Found the solution
        result = list(perm)
        print(result)
        break