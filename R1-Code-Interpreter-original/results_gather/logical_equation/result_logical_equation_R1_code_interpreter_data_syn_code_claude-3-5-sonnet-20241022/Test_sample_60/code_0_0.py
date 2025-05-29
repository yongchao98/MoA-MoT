from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if A - G != -13: return False
    if B + C != 15: return False
    if abs(H - 4.8*B) > 0.1: return False  # Using small tolerance for floating point
    if I - B != -3: return False
    if L - J != -5: return False
    if B <= A: return False
    if abs(I - 2.0*F) > 0.1: return False  # Using small tolerance for floating point
    if C - H != -14: return False
    if B - I != 3: return False
    if A - F != 2: return False
    
    return True

# Numbers to assign
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break