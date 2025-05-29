def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if A - D != -77: return False
    if H <= A: return False
    if I != 3 * E: return False
    if I <= B: return False
    if E - H != -27: return False
    if G - A != -2: return False
    if H != 4 * E: return False
    
    return True

def find_solution():
    numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
    
    # We already know some values
    E = 9
    H = 36
    I = 27
    
    # Try all remaining combinations
    for A in numbers:
        if A >= H: continue  # H > A must be satisfied
        
        D = A + 77
        if D not in numbers: continue
        
        G = A - 2
        if G not in numbers: continue
        
        remaining = [n for n in numbers if n not in [A, D, E, G, H, I]]
        
        for B in remaining:
            if I <= B: continue  # I > B must be satisfied
            
            for C in [n for n in remaining if n != B]:
                for F in [n for n in remaining if n not in [B, C]]:
                    perm = [A, B, C, D, E, F, G, H, I]
                    if len(set(perm)) == 9 and check_conditions(perm):
                        return perm

solution = find_solution()
print(solution)