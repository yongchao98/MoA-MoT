def find_solution():
    # Known values from direct equations:
    # From I - J = 51 and I + J = 141:
    # 2I = 192 -> I = 96
    # 2J = 90 -> J = 45
    I = 96
    J = 45
    
    # From H + I = 103:
    # H + 96 = 103 -> H = 7
    H = 7
    
    # From E = 4.0H:
    # E = 4.0 * 7 = 28
    E = 28
    
    # From E - D = 25:
    # 28 - D = 25 -> D = 3
    D = 3
    
    # From E = 2.8M:
    # 28 = 2.8M -> M = 10
    M = 10
    
    # From G - E = -13:
    # G - 28 = -13 -> G = 15
    G = 15
    
    # From G = 3.0A:
    # 15 = 3.0A -> A = 5
    A = 5
    
    # From C + I = 120:
    # C + 96 = 120 -> C = 24
    C = 24
    
    # Remaining numbers: 2, 16, 36, 50
    # Must be assigned to B, F, K, L
    
    # B = 2 (smallest remaining)
    B = 2
    
    # F = 16 (next smallest)
    F = 16
    
    # K = 36, L = 50 (remaining two in ascending order)
    K = 36
    L = 50
    
    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
    print(result)

find_solution()