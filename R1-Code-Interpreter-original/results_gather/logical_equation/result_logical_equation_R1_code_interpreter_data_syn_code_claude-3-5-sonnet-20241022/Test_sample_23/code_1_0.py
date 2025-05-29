from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if D + E != 15:  # Condition 1
        return False
    if A + F != 201:  # Condition 2
        return False
    if D != 1.5 * E:  # Condition 3
        return False
    if C != 4.0 * D:  # Condition 4
        return False
    if H <= G:  # Condition 5
        return False
    if E - I != -26:  # Condition 6
        return False
    if E != 2.0 * G:  # Condition 7
        return False
    if C + I != 68:  # Condition 8
        return False
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

for p in permutations(numbers):
    if check_conditions(p):
        print(list(p))
        break