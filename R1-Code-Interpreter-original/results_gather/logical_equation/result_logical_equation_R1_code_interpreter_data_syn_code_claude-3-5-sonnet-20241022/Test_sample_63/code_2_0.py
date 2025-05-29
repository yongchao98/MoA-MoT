from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    if A - E != 4: return False
    if E - A != -4: return False
    if I != 4.0 * J: return False
    if D - B != -14: return False
    if H + M != 26: return False
    if L != int(3.2 * A): return False
    if L != int(1.6 * F): return False
    if F + L != 26: return False  # This is crucial
    if H - I != -4: return False
    if D != int(1.5 * H): return False
    
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break