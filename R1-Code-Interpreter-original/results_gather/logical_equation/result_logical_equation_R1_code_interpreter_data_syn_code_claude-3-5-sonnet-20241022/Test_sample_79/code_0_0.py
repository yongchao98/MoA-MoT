from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if F + J != 15: return False  # 1
    if F - G != -91: return False  # 2
    if abs(D - 3.5*K) > 0.01: return False  # 3
    if E + G != 141: return False  # 4
    if A - K != 13: return False  # 5
    if F + M != 29: return False  # 6
    if abs(D - 1.4*F) > 0.01: return False  # 7
    if abs(L - 4.0*D) > 0.01: return False  # 8
    if abs(C - 3.6*J) > 0.01: return False  # 9
    if A + M != 39: return False  # 10
    if abs(C - 2.4*A) > 0.01: return False  # 11
    
    return True

# Numbers to assign
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all possible permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break