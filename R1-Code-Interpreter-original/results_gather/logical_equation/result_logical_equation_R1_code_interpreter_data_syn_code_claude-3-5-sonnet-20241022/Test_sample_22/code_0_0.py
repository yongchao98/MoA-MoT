from itertools import permutations

def check_conditions(nums):
    A, B, C, D, E, F, G, H, I = nums
    
    # Check all conditions
    if B != 2.0 * G:  # Condition 1
        return False
    if I != 2.5 * F:  # Condition 2
        return False
    if B + I != 86:   # Condition 3
        return False
    if A + G != 12:   # Condition 4
        return False
    if I <= F:        # Condition 5
        return False
    if F + I != 112:  # Condition 6
        return False
    if C + H != 140:  # Condition 7
        return False
    if A - H != -117: # Condition 8
        return False
    return True

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break