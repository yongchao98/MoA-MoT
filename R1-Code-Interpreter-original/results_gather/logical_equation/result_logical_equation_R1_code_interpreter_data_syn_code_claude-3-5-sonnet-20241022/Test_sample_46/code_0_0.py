from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K = [values[letters.index(x)] for x in 'ABCDEFGHIJK']
    
    # Check all constraints
    if A + F != 196: return False
    if G - D != -6: return False
    if G - A != -72: return False
    if E != 3.0 * A: return False
    if B <= I: return False
    if K - A != 75: return False
    if H - F != -41: return False
    if B - G != 36: return False
    if A - H != -5: return False
    if E - A != 150: return False
    if H - J != 64: return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
letters = list('ABCDEFGHIJK')

for perm in permutations(numbers):
    if check_constraints(perm, letters):
        result = list(perm)
        print(result)
        break