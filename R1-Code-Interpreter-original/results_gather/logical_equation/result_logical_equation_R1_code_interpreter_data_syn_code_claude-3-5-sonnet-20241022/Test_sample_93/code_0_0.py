from itertools import permutations

def check_constraints(values, mapping):
    A,B,C,D,E,F,G,H,I,J,K,L,M = values
    
    # Check all constraints
    if C + E != 141: return False  # 1
    if J + L != 33: return False   # 2
    if abs(K - 2.4*H) > 0.1: return False  # 3 (using small tolerance)
    if B - J != 10: return False   # 4
    if E - C != 51: return False   # 5
    if C <= F: return False        # 6
    if abs(M - 3.2*J) > 0.1: return False  # 7 (using small tolerance)
    if E + K != 120: return False  # 8
    if J - H != -5: return False   # 9
    if K - B != 9: return False    # 10
    if B - M != -1: return False   # 11
    if C <= M: return False        # 12
    
    return True

# Possible values
values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for perm in permutations(values):
    if check_constraints(perm, values):
        print(list(perm))
        break