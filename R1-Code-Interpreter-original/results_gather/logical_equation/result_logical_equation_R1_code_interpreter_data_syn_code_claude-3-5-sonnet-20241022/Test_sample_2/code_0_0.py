from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if G - E != -3:  # Condition 1
        return False
    if D - C != -41:  # Condition 2
        return False
    if C <= B:  # Condition 3
        return False
    if F != 4.0 * H:  # Condition 4
        return False
    if I - G != 2:  # Condition 5
        return False
    if A + E != 31:  # Condition 6
        return False
    if A - B != 11:  # Condition 7
        return False
    
    return True

numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break