from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I = [values[letters.index(c)] for c in 'ABCDEFGHI']
    
    # Check all conditions with exact equality
    if abs(G - 1.5 * H) > 0.0001:  # G = 1.5H
        return False
    if abs((E - D) - (-41)) > 0.0001:  # E - D = -41
        return False
    if abs(F - 1.5 * E) > 0.0001:  # F = 1.5E
        return False
    if abs((D + G) - 130) > 0.0001:  # D + G = 130
        return False
    if abs((F + H) - 126) > 0.0001:  # F + H = 126
        return False
    if abs((F + I) - 123) > 0.0001:  # F + I = 123
        return False
    if abs(E - 2.5 * B) > 0.0001:  # E = 2.5B
        return False
    if abs((D - E) - 41) > 0.0001:  # D - E = 41
        return False
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = list('ABCDEFGHI')

for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(perm)