from itertools import permutations

def check_solution(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if A != 3.0 * E:  # A = 3.0E
        return False
    if D + F != 81:   # D + F = 81
        return False
    if A + D != 28:   # A + D = 28
        return False
    if C + F != 96:   # C + F = 96
        return False
    if A + E != 36:   # A + E = 36
        return False
    if B + E != 12:   # B + E = 12
        return False
    if B - F != -77:  # B - F = -77
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_solution(perm):
        print(list(perm))
        break