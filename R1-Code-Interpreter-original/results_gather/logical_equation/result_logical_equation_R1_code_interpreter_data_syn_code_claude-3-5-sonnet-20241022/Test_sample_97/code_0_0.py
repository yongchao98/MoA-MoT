from itertools import permutations

def check_conditions(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all conditions
    if A - L != -5:  # Condition 1
        return False
    if abs(I - 3.2 * C) > 0.1:  # Condition 2 (with small tolerance)
        return False
    if J + M != 60:  # Condition 3
        return False
    if L + M != 31:  # Condition 4
        return False
    if abs(C - 2.5 * A) > 0.1:  # Condition 5
        return False
    if abs(L - 3.5 * A) > 0.1:  # Condition 6
        return False
    if abs(D - 4.0 * M) > 0.1:  # Condition 7
        return False
    if abs(B - 4.5 * G) > 0.1:  # Condition 8
        return False
    if I + M != 40:  # Condition 9
        return False
    if A + D != 98:  # Condition 10
        return False
    if abs(F - 1.5 * G) > 0.1:  # Condition 11
        return False
    if D + G != 106:  # Condition 12
        return False
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break