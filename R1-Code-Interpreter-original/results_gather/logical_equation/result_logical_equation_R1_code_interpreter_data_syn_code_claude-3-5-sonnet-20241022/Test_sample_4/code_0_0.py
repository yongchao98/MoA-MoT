from itertools import permutations

def check_solution(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all constraints
    if C - E != -5:  # 1
        return False
    if D - F != 79:  # 2
        return False
    if G != 4.0 * E:  # 3
        return False
    if H != 4.0 * C:  # 4
        return False
    if E - G != -27:  # 5
        return False
    if E + H != 25:  # 6
        return False
    if A + H != 137:  # 7
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_solution(perm):
        print(list(perm))
        break