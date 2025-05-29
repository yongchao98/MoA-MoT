from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if not (D > A): return False
    if not (D - A == 66): return False
    if not (F - G == -145): return False
    if not (B - G == -209): return False
    if not (abs(E - 3.0 * B) < 0.1): return False  # Using small tolerance for float comparison
    if not (I - A == 112): return False
    if not (D > C): return False
    if not (A + H == 159): return False
    if not (D + H == 225): return False
    if not (abs(E - 2.4 * J) < 0.1): return False  # Using small tolerance for float comparison
    if not (A - H == -141): return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break