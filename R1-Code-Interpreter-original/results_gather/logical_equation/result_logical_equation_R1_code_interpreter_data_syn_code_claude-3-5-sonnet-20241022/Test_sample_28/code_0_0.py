from itertools import permutations

def check_conditions(values, positions):
    A, B, C, D, E, F, G, H, I = [values[positions[c]] for c in 'ABCDEFGHI']
    
    # Check all conditions
    if G != 1.5 * C:  # G = 1.5C
        return False
    if E != 1.8 * D:  # E = 1.8D
        return False
    if D + I != 141:  # D + I = 141
        return False
    if B - H != -88:  # B - H = -88
        return False
    if H + I != 241:  # H + I = 241
        return False
    if H - G != 111:  # H - G = 111
        return False
    if H - B != 88:   # H - B = 88
        return False
    if A + I != 201:  # A + I = 201
        return False
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = 'ABCDEFGHI'

# Try all possible permutations
for perm in permutations(range(9)):
    positions = {letters[i]: perm[i] for i in range(9)}
    if check_conditions(numbers, positions):
        result = [numbers[positions[c]] for c in 'ABCDEFGHI']
        print(result)
        break