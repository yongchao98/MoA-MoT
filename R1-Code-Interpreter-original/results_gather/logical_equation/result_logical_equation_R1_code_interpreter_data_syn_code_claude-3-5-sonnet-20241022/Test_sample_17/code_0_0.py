from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all conditions
    if abs(A - 1.6 * G) > 0.1:  # Allow small floating point error
        return False
    if E + H != 129:
        return False
    if abs(H - 1.5 * B) > 0.1:  # Allow small floating point error
        return False
    if E + I != 130:
        return False
    if A + B != 112:
        return False
    if D + G != 23:
        return False
    if abs(B - 2.5 * A) > 0.1:  # Allow small floating point error
        return False
    if abs(B - 4.0 * G) > 0.1:  # Allow small floating point error
        return False
    
    return True

numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break