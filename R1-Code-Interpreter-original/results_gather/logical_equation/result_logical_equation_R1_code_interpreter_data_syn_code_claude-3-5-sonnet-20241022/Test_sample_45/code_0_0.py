from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if F != 2.4 * D:  # F = 2.4D
        return False
    if G != 1.5 * K:  # G = 1.5K
        return False
    if G + K != 375:  # G + K = 375
        return False
    if F - G != -177:  # F - G = -177
        return False
    if A + D != 29:  # A + D = 29
        return False
    if I - D != 19:  # I - D = 19
        return False
    if A != 3.0 * H:  # A = 3.0H
        return False
    if F != 3.0 * B:  # F = 3.0B
        return False
    if F + K != 198:  # F + K = 198
        return False
    if G <= C:  # G > C
        return False
    if J - F != 73:  # J - F = 73
        return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break