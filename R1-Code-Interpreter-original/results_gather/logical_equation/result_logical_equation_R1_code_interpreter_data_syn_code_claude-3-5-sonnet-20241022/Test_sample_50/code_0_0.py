from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    if B + F != 346: return False
    if A + E != 155: return False
    if F + G != 228: return False
    if A + H != 128: return False
    if abs(A - 4.0 * D) > 0.01: return False  # Using small epsilon for float comparison
    if I - K != -23: return False
    if abs(J - 2.0 * E) > 0.01: return False
    if abs(F - 3.0 * E) > 0.01: return False
    if B - J != -29: return False
    if A + K != 119: return False
    if abs(H - 3.0 * I) > 0.01: return False
    
    return True

# Numbers to assign
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break