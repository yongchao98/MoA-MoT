from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if not B > F:  # Condition 1
        return False
    if D + E != 107:  # Condition 2
        return False
    if not D > F:  # Condition 3
        return False
    if D - F != 77:  # Condition 4
        return False
    if not C > B:  # Condition 5
        return False
    if C + I != 52:  # Condition 6
        return False
    if abs(B - 3.0 * F) > 0.0001:  # Condition 7, using small epsilon for float comparison
        return False
    
    return True

numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
valid_solutions = []

for perm in permutations(numbers):
    if check_conditions(perm):
        valid_solutions.append(perm)

if valid_solutions:
    print(list(valid_solutions[0]))
else:
    print("No solution found")