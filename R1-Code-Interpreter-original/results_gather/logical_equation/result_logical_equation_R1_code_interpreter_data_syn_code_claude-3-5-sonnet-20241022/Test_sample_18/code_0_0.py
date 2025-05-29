from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all conditions
    if A + C != 127:  # Condition 1
        return False
    if D - E != -77:  # Condition 2
        return False
    if A - B != 89:  # Condition 3
        return False
    if abs(B - 1.6 * F) > 0.1:  # Condition 4 (using small tolerance for float comparison)
        return False
    if E <= I:  # Condition 5
        return False
    if abs(E - 4.0 * F) > 0.1:  # Condition 6
        return False
    if abs(I - 4.0 * H) > 0.1:  # Condition 7
        return False
    if abs(H - 1.5 * C) > 0.1:  # Condition 8
        return False
    
    return True

# Numbers to assign
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break