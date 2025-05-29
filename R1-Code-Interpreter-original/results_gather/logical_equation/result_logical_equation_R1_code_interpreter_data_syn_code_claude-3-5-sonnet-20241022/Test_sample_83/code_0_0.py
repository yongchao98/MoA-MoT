def check_equations(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all equations
    if E - D != 25: return False
    if H + I != 103: return False
    if abs(E - 2.8*M) > 0.1: return False  # Using small tolerance for float comparison
    if C + I != 120: return False
    if A + I != 101: return False
    if I - J != 51: return False
    if I + J != 141: return False
    if H - I != -89: return False
    if abs(E - 4.0*H) > 0.1: return False
    if abs(G - 3.0*A) > 0.1: return False
    if G - E != -13: return False
    
    return True

def solve():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    from itertools import permutations
    
    # From I - J = 51 and I + J = 141, we can deduce:
    # I = 96 and J = 45
    
    for perm in permutations(numbers, 13):
        if perm[8] != 96 or perm[9] != 45:  # I and J must be 96 and 45
            continue
            
        if check_equations(perm):
            return list(perm)
    
    return None

result = solve()
if result:
    print(result)