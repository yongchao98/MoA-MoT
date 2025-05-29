from itertools import permutations

def check_constraints(values, numbers):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all constraints
    if C + I != 241:  # Constraint 1
        return False
    if E - H != -44:  # Constraint 2
        return False
    if I != 1.5 * H:  # Constraint 3
        return False
    if D + E != 42:   # Constraint 4
        return False
    if D + H != 86:   # Constraint 5
        return False
    if F + G != 35:   # Constraint 6
        return False
    if D != 2.0 * F:  # Constraint 7
        return False
    if A + I != 140:  # Constraint 8
        return False
    
    return True

# Numbers that can be assigned
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm, numbers):
        print(list(perm))
        break