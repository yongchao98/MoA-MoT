from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Equations and inequalities
def check_solution(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all conditions
    if not (L - I == -23):
        return False
    if not (I + M == 44):
        return False
    if not (D + K == 98):
        return False
    if not (F + I == 73):
        return False
    if not (F > J):
        return False
    if not (B - C == 26):
        return False
    if not (B == 2.4 * A):
        return False
    if not (C + F == 55):
        return False
    if not (F - J == 38):
        return False
    if not (H - F == 5):
        return False
    if not (J == 1.4 * L):
        return False
    
    return True

# Try all permutations of the numbers
for perm in permutations(numbers):
    if check_solution(perm):
        print(list(perm))
        break