from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if I != 4.0 * C:  # I = 4.0C
        return False
    if F + I != 83:  # F + I = 83
        return False
    if K - E != 82:  # K - E = 82
        return False
    if E + F != 42:  # E + F = 42
        return False
    if A + H != 225:  # A + H = 225
        return False
    if B + K != 130:  # B + K = 130
        return False
    if A <= B:  # A > B
        return False
    if J != 3.0 * H:  # J = 3.0H
        return False
    if B != 3.0 * F:  # B = 3.0F
        return False
    if C + F != 23:  # C + F = 23
        return False
    if A <= K:  # A > K
        return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break