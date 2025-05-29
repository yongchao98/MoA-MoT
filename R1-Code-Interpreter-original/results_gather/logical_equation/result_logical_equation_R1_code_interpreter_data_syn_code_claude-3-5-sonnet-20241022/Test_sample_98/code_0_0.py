from itertools import permutations

def check_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if abs(A - 2.8*K) > 0.1:  # Allow small floating point error
        return False
    if A + F != 30:
        return False
    if C + D != 31:
        return False
    if G - L != 45:
        return False
    if L - C != -19:
        return False
    if abs(I - 1.5*K) > 0.1:  # Allow small floating point error
        return False
    if A <= D:
        return False
    if J + K != 81:  # 36 + 45 = 81
        return False
    if abs(B - 3.6*K) > 0.1:  # Allow small floating point error
        return False
    if abs(C - 1.6*I) > 0.1:  # Allow small floating point error
        return False
    if B + G != 86:
        return False
    if K - J != 9:  # 45 - 36 = 9
        return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break