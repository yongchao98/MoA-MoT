from itertools import permutations

def check_constraints(p):
    A, B, C, D, E, F, G, H, I = p
    
    # Check all constraints
    if C != 3.0 * B:  # C = 3.0B
        return False
    if A != 4.0 * C:  # A = 4.0C
    	return False
    if F <= A:  # F > A
        return False
    if D != 4.0 * H:  # D = 4.0H
        return False
    if E <= B:  # E > B
        return False
    if A - H != 32:  # A - H = 32
        return False
    if F + H != 84:  # F + H = 84
        return False
    
    return True

numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break