from itertools import permutations

def check_constraints(p):
    A, B, C, D, E, F, G, H, I = p
    
    # Check all constraints
    if D - E != 24:  # Constraint 1
        return False
    if A - F != 35:  # Constraint 2
        return False
    if D != 3.0 * B:  # Constraint 3
        return False
    if I != 4.0 * F:  # Constraint 4
        return False
    if B != 3.0 * E:  # Constraint 5
        return False
    if E + H != 124:  # Constraint 6
        return False
    if C <= E:  # Constraint 7
        return False
    
    return True

# Numbers that can be assigned
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Try all permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break