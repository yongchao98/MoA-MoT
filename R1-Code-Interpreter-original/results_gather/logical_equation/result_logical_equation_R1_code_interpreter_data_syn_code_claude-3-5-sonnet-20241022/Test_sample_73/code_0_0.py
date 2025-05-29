from itertools import permutations

def check_conditions(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all conditions
    if M != 3.0 * L:  # M = 3.0L
        return False
    if F != 2.5 * H:  # F = 2.5H
        return False
    if E - L != 9:  # E - L = 9
        return False
    if A + I != 60:  # A + I = 60
        return False
    if E + J != 17:  # E + J = 17
        return False
    if G - I != 21:  # G - I = 21
        return False
    if G - A != -9:  # G - A = -9
        return False
    if E + G != 46:  # E + G = 46
        return False
    if F + M != 8:  # F + M = 8
        return False
    if H + I != 17:  # H + I = 17
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break