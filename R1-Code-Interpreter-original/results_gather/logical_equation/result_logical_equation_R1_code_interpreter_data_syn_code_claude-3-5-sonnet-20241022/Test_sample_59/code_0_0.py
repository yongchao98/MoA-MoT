from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if not (K > A): return False
    if abs(H - 2.4*A) > 0.1: return False  # Using small tolerance for float comparison
    if abs(H - 3.0*F) > 0.1: return False
    if B - I != -216: return False
    if C - J != -72: return False
    if abs(I - 1.5*G) > 0.1: return False
    if abs(I - 3.0*J) > 0.1: return False
    if B - H != -39: return False
    if D - I != -104: return False
    if J - F != 59: return False
    if B - F != -7: return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break