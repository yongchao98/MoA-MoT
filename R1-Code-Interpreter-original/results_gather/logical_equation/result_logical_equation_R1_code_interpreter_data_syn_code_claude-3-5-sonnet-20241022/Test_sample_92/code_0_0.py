from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if not abs(A - 3.0 * L) < 0.1:  # Condition 1: A = 3.0L
        return False
    if A + M != 39:  # Condition 2: A + M = 39
        return False
    if L + M != 29:  # Condition 3: L + M = 29
        return False
    if not abs(B - 3.2 * L) < 0.1:  # Condition 4: B = 3.2L
        return False
    if not abs(A - 1.5 * I) < 0.1:  # Condition 5: A = 1.5I
        return False
    if C - B != 80:  # Condition 6: C - B = 80
        return False
    if E + K != 43:  # Condition 7: E + K = 43
        return False
    if C <= D:  # Condition 8: C > D
        return False
    if not abs(E - 2.4 * A) < 0.1:  # Condition 9: E = 2.4A
        return False
    if C + D != 124:  # Condition 10: C + D = 124
        return False
    if H - A != -12:  # Condition 11: H - A = -12
        return False
    if F + G != 47:  # Condition 12: F + G = 47
        return False
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break