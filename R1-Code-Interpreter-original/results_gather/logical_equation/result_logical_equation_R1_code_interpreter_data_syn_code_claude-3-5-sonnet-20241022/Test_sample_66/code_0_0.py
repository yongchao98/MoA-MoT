from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if D != 3.0 * I:  # D = 3.0I
        return False
    if C + H != 20:  # C + H = 20
        return False
    if A <= F:  # A > F
        return False
    if F + G != 38:  # F + G = 38
        return False
    if L - G != -26:  # L - G = -26
        return False
    if C + F != 17:  # C + F = 17
        return False
    if E + L != 60:  # E + L = 60
        return False
    if M != 4.5 * L:  # M = 4.5L
        return False
    if A - K != -17:  # A - K = -17
        return False
    if D + H != 8:  # D + H = 8
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break