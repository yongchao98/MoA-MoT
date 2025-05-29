from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if D + J != 252:  # Constraint 1
        return False
    if abs(D - 3.0 * K) > 0.01:  # Constraint 2 (using small epsilon for float comparison)
        return False
    if H - A != 44:  # Constraint 3
        return False
    if J <= E:  # Constraint 4
        return False
    if abs(C - 2.0 * K) > 0.01:  # Constraint 5
        return False
    if D <= A:  # Constraint 6
        return False
    if A + F != 52:  # Constraint 7
        return False
    if D + E != 228:  # Constraint 8
        return False
    if abs(B - 3.0 * F) > 0.01:  # Constraint 9
        return False
    if E + I != 124:  # Constraint 10
        return False
    
    return True

# Numbers to assign
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break