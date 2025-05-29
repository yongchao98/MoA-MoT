from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if D + E != 52:  # Condition 1
        return False
    if A + F != 201:  # Condition 2
        return False
    if D * 2 != E * 3:  # Condition 3 (D = 1.5E)
        return False
    if C != 4 * D:  # Condition 4
        return False
    if H <= G:  # Condition 5
        return False
    if E - I != -89:  # Condition 6
        return False
    if E != 2 * G:  # Condition 7
        return False
    if C + I != 141:  # Condition 8
        return False
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
found = False

for p in permutations(numbers):
    if check_conditions(p):
        print(list(p))
        found = True
        break

if not found:
    print("No solution found")