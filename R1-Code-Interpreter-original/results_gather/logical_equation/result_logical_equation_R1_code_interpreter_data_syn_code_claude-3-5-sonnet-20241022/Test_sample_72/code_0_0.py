from itertools import permutations

def check_constraints(perm):
    # Map letters to values
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all constraints
    if F + J != 86: return False  # 1
    if B != 4.5 * A: return False  # 2
    if A + K != 15: return False  # 3
    if M != 3.0 * L: return False  # 4
    if H != 1.5 * A: return False  # 5
    if J <= D: return False  # 6
    if A - M != 7: return False  # 7
    if B + E != 52: return False  # 8
    if A <= G: return False  # 9
    if H <= L: return False  # 10
    
    return True

# Numbers to assign
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        print(list(perm))
        break