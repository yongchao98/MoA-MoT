from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if A - C != 26:  # Condition 1
        return False
    if M <= L:  # Condition 2
        return False
    if abs(C - 2.0 * D) > 0.01:  # Condition 3
        return False
    if E + I != 57:  # Condition 4
        return False
    if abs(L - 1.6 * C) > 0.01:  # Condition 5
        return False
    if G - K != -12:  # Condition 6
        return False
    if abs(K - 1.5 * C) > 0.01:  # Condition 7
        return False
    if abs(A - 3.6 * C) > 0.01:  # Condition 8
        return False
    if abs(M - 4.8 * D) > 0.01:  # Condition 9
        return False
    if C + F != 106:  # Condition 10
        return False
    if H - G != 42:  # Condition 11
        return False
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break