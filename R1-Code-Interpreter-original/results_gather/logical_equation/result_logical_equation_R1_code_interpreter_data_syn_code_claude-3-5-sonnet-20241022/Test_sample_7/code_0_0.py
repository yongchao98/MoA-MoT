from itertools import permutations

def check_conditions(values, numbers):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all conditions
    if D + H != 13:
        return False
    if E - C != 15:
        return False
    if C + H != 5:
        return False
    if H != 4.0 * C:
        return False
    if H - B != -117:
        return False
    if F != 3.0 * C:
        return False
    if F + H != 7:
        return False
    
    return True

# The possible numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm, numbers):
        print(list(perm))