from itertools import permutations

def check_equations(values, letters):
    A, B, C, D, E, F, G, H, I = [values[letters.index(c)] for c in 'ABCDEFGHI']
    
    # Check all equations
    if H != 1.6 * G:  # H = 1.6G
        return False
    if F + I != 9:    # F + I = 9
        return False
    if D != 2.5 * H:  # D = 2.5H
        return False
    if B + D != 201:  # B + D = 201
        return False
    if I - D != -77:  # I - D = -77
        return False
    if F != 2.0 * I:  # F = 2.0I
        return False
    if C + D != 89:   # C + D = 89
        return False
    if E + I != 39:   # E + I = 39
        return False
    
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = list('ABCDEFGHI')

# Try all permutations
for perm in permutations(numbers):
    if check_equations(perm, letters):
        print(list(perm))
        break