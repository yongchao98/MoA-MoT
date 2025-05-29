from itertools import permutations

def check_constraints(p):
    A, B, C, D, E, F, G, H, I, J, K = p
    
    # Check all constraints
    if J + K != 271: return False  # 1
    if G <= I: return False        # 2
    if E + J != 141: return False  # 3
    if A + C != 300: return False  # 4
    if G + K != 230: return False  # 5
    if B - K != -141: return False # 6
    if D - C != -209: return False # 7
    if A + I != 78: return False   # 8
    if F - D != 23: return False   # 9
    if E - H != -28: return False  # 10
    if H + K != 198: return False  # 11
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all permutations
for p in permutations(numbers):
    if check_constraints(p):
        print(list(p))
        break