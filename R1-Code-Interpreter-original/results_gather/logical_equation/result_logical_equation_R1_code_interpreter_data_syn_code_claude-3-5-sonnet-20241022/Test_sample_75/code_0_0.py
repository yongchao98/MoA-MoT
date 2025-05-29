from itertools import permutations

def check_conditions(values, perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if L - I != -23:  # Condition 1
        return False
    if I + M != 44:   # Condition 2
        return False
    if D + K != 98:   # Condition 3
        return False
    if F + I != 73:   # Condition 4
        return False
    if F <= J:        # Condition 5
        return False
    if B - C != 26:   # Condition 6
        return False
    if abs(B - 2.4*A) > 0.0001:  # Condition 7 (using small epsilon for float comparison)
        return False
    if C + F != 55:   # Condition 8
        return False
    if F - J != 38:   # Condition 9
        return False
    if H - F != 5:    # Condition 10
        return False
    if abs(J - 1.4*L) > 0.0001:  # Condition 11
        return False
    return True

# The possible values
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(values):
    if check_conditions(values, perm):
        print(list(perm))
        break