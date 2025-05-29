def check_equations(F, L):
    # Given F and L, check if we can satisfy all equations
    numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
    used = {F, L}
    
    # B - F = -21
    B = F - 21
    if B not in numbers or B in used:
        return None
    used.add(B)
    
    # K - L = -42
    K = L - 42
    if K not in numbers or K in used:
        return None
    used.add(K)
    
    # F - A = -22
    A = F + 22
    if A not in numbers or A in used:
        return None
    used.add(A)
    
    # B + C = 43
    C = 43 - B
    if C not in numbers or C in used:
        return None
    used.add(C)
    
    # D + K = 4
    D = 4 - K
    if D not in numbers or D in used:
        return None
    used.add(D)
    
    # H + L = 60
    H = 60 - L
    if H not in numbers or H in used:
        return None
    used.add(H)
    
    # H + M = 17
    M = 17 - H
    if M not in numbers or M in used:
        return None
    used.add(M)
    
    # C = 2.4H
    if abs(C - 2.4*H) > 0.01:
        return None
    
    # G = 2.5M
    G = int(2.5 * M)
    if G not in numbers or G in used:
        return None
    used.add(G)
    
    # Find remaining unused numbers for E, I, J
    remaining = [x for x in numbers if x not in used]
    if len(remaining) != 3:
        return None
    
    E, I, J = remaining
    
    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
    return result

# Try all possible F, L combinations where F + L = 73
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
for F in numbers:
    for L in numbers:
        if F != L and F + L == 73:
            result = check_equations(F, L)
            if result is not None:
                print(result)
                exit()