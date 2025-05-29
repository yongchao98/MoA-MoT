from itertools import permutations

def check_conditions(values, mapping):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if not E > C:  # 1
        return False
    if abs(G - 4.8 * M) > 0.1:  # 2 (using small tolerance for float comparison)
        return False
    if L - C != 12:  # 3
        return False
    if D + H != 66:  # 4
        return False
    if L - D != -35:  # 5
        return False
    if abs(E - 3.6 * B) > 0.1:  # 6
        return False
    if B - D != -40:  # 7
        return False
    if H - G != -8:  # 8
        return False
    if abs(C - 1.5 * K) > 0.1:  # 9
        return False
    if abs(H - 1.6 * B) > 0.1:  # 10
        return False
    if J - E != 9:  # 11
        return False
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        break