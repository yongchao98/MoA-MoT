from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I = [values[letters.index(c)] for c in 'ABCDEFGHI']
    
    # Check all conditions
    if G != 1.5 * H:  # G = 1.5H
        return False
    if E - D != -41:  # E - D = -41
        return False
    if F != 1.5 * E:  # F = 1.5E
        return False
    if D + G != 130:  # D + G = 130
        return False
    if F + H != 126:  # F + H = 126
        return False
    if F + I != 123:  # F + I = 123
        return False
    if E != 2.5 * B:  # E = 2.5B
        return False
    if D - E != 41:   # D - E = 41
        return False
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = list('ABCDEFGHI')
found = False

for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        found = True
        break

if not found:
    print("No solution found")