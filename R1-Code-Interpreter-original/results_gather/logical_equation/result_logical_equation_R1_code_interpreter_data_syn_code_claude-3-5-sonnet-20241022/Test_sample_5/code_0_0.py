from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if G != 4.0 * A:  # Condition 1
        return False
    if G <= A:  # Condition 2
        return False
    if C != 4.0 * D:  # Condition 3
        return False
    if B + E != 148:  # Condition 4
        return False
    if A - E != -117:  # Condition 5
        return False
    if D - F != -71:  # Condition 6
        return False
    if C + H != 39:  # Condition 7
        return False
    
    return True

numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
solution = None

for perm in permutations(numbers):
    if check_conditions(perm):
        solution = perm
        break

if solution:
    print(list(solution))