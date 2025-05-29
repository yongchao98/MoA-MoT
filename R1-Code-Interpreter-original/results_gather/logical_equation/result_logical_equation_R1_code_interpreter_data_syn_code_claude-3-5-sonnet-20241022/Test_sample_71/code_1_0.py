from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Quick checks first to fail fast
    if B - F != -21: return False
    if F + L != 73: return False  # F and L sum must be 73
    
    # Pre-filter based on known relationships
    if F - A != -22: return False
    if K - L != -42: return False
    if B + C != 43: return False
    if H + L != 60: return False
    if H + M != 17: return False
    
    # Float comparisons
    if abs(C - 2.4*H) > 0.01: return False
    if D + K != 4: return False
    if abs(G - 2.5*M) > 0.01: return False
    
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Pre-filter possible F and L combinations (F + L = 73)
possible_F_L = []
for f in numbers:
    for l in numbers:
        if f + l == 73:
            possible_F_L.append((f, l))

# For each valid F,L pair, try to find other values
for F, L in possible_F_L:
    remaining_numbers = [n for n in numbers if n not in [F, L]]
    
    # Find possible B values (B - F = -21)
    B = F - 21
    if B not in remaining_numbers:
        continue
        
    remaining_numbers.remove(B)
    
    # Find possible K values (K - L = -42)
    K = L - 42
    if K not in remaining_numbers:
        continue
        
    remaining_numbers.remove(K)
    
    # Try permutations of remaining numbers for other positions
    positions = [0, 2, 3, 4, 6, 7, 8, 9]  # Indices for A,C,D,E,G,H,I,J,M
    for perm in permutations(remaining_numbers, len(positions)):
        values = [None] * 13
        for pos, val in zip(positions, perm):
            values[pos] = val
        values[1] = B  # B
        values[5] = F  # F
        values[10] = K  # K
        values[11] = L  # L
        
        if check_conditions(values):
            print(values)
            exit()