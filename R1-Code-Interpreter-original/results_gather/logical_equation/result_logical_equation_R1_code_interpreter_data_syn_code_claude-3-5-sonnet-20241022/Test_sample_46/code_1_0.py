def find_solution():
    # Known values from direct equations
    # E = 3A and E - A = 150 means:
    # 3A - A = 150
    # 2A = 150
    # A = 75
    A = 75
    E = 225  # E = 3A = 3*75 = 225
    
    # K - A = 75, so:
    K = 150  # K = A + 75 = 75 + 75 = 150
    
    # G - A = -72, so:
    G = 3    # G = A - 72 = 75 - 72 = 3
    
    # A + F = 196, so:
    F = 121  # F = 196 - A = 196 - 75 = 121
    
    # H - F = -41, so:
    H = 80   # H = F - 41 = 121 - 41 = 80
    
    # G - D = -6, so:
    D = 9    # D = G + 6 = 3 + 6 = 9
    
    # B - G = 36, so:
    B = 39   # B = G + 36 = 3 + 36 = 39
    
    # H - J = 64, so:
    J = 16   # J = H - 64 = 80 - 64 = 16
    
    # Only remaining number for I must be 20
    # This satisfies B > I since 39 > 20
    I = 20
    
    # C must be 16 (the only remaining number)
    C = 16
    
    result = [A, B, C, D, E, F, G, H, I, J, K]
    print(result)

find_solution()